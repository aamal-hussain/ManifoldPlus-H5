import h5py
import pyvista as pv
import numpy as np
import os
import shutil
import signal
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from tempfile import TemporaryDirectory
import subprocess
import time
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    filename='manifold_processing_test.log'
)

DATASET_SAVE_DIR = "/mnt/storage01/workspace/research/lethe/02_training_datasets/shapenet_50k_manifold_core/test/"
FAILED_SAVE_DIR = "./shapenet_50k_failed/"
RAW_DATA_DIR = "/mnt/storage01/workspace/research/lethe/02_training_datasets/shapenet_core_v2/test/"
MAX_RETRIES = 3
TIMEOUT_SECONDS = 300

def outward_point_normal_check(mesh: pv.PolyData):
    """Simple check to ensure normals are oriented correctly (outward).
    Looks at the point with the largest value in Z, and checks the Z component of the normal
    is positive.
    Args:
        mesh (pv.PolyData): mesh with point data
    Raises:
        AttributeError: if the normals are not oriented correctly
    """
    max_z_index = mesh.points[:, 2].argmax()
    if mesh.point_data["Normals"][max_z_index][2] < 0:
        raise AttributeError("Doubtful the normals are oriented correctly!")

def outward_cell_normal_check(mesh: pv.PolyData):
    """Simple check to ensure normals are oriented correctly (outward).
    Finds all the faces containing the 'highest' point (ie. the point with the largest Z value).
    Of those faces, finds the uppermost one (ie. has the largest Z component in its normal).
    If the Z component of the normal of the uppermost face is negative, then raises an exception.
    Args:
        mesh (pv.PolyData): mesh with cell data
    Raises:
        AttributeError: if the normals are not oriented correctly
    """
    index_of_max_z_point = mesh.points[:, 2].argmax()
    mesh_containg_max_z_point = mesh.extract_points(index_of_max_z_point)
    normals = mesh_containg_max_z_point.cell_data["Normals"]
    index_of_uppermost_face = np.argmax(np.abs(normals[:, 2]))
    if normals[index_of_uppermost_face, 2] < 0:
        raise AttributeError("Doubtful the normals are oriented correctly!")

def run_manifold_command(input_path, output_path):
    """Run manifold command with timeout and proper error handling."""
    try:
        process = subprocess.Popen(
            ["../build/manifold", "--source", input_path, "--target", output_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            preexec_fn=os.setsid 
        )
        
        stdout, stderr = process.communicate(timeout=TIMEOUT_SECONDS)
        
        if process.returncode != 0:
            raise subprocess.CalledProcessError(
                process.returncode, 
                ["../build/manifold"], 
                output=stdout, 
                stderr=stderr
            )
            
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        process.kill()
        process.wait()
        raise TimeoutError(f"Process timed out after {TIMEOUT_SECONDS} seconds")
        
    except Exception as e:
        try:
            os.killpg(os.getpgid(process.pid), signal.SIGTERM)
            process.kill()
            process.wait()
        except:
            pass
        raise e

def make_manifold(file):
    input_path = RAW_DATA_DIR + file
    retry_count = 0
    last_exception = None
    
    while retry_count < MAX_RETRIES:
        try:
            with TemporaryDirectory() as tmpdir:
                output_path = f"{tmpdir}/{file}"
                run_manifold_command(input_path, output_path)
                
                if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
                    raise ValueError("Output file is empty or doesn't exist")
                
                with h5py.File(output_path, "r") as f:
                    verts = f["mesh.verts"][:]
                    faces = f["mesh.faces"][:]
                    verts_normals = f["mesh.verts_normals"][:]
                    faces_normals = f["mesh.faces_normals"][:]
                    
                    mesh = pv.PolyData.from_regular_faces(verts, faces)
                    mesh.point_data["Normals"] = verts_normals
                    mesh.cell_data["Normals"] = faces_normals
                    mesh = mesh.compute_cell_sizes(length=False, area=True, volume=False)
                
                    outward_cell_normal_check(mesh)
                    outward_point_normal_check(mesh)
                
                if mesh.is_manifold:
                    shutil.move(output_path, DATASET_SAVE_DIR + file)
                else:
                    shutil.move(output_path, FAILED_SAVE_DIR + file)
                
                return mesh.is_manifold
                
        except (subprocess.CalledProcessError, TimeoutError, ValueError, AttributeError) as e:
            last_exception = e
            retry_count += 1
            logging.warning(f"Attempt {retry_count} failed for {file}: {str(e)}")
            time.sleep(1) 
            continue
            
        except Exception as e:
            logging.error(f"Unrecoverable error in {file}: {str(e)}")
            raise e
            
    logging.error(f"All {MAX_RETRIES} attempts failed for {file}. Last error: {str(last_exception)}")
    return False

def manifold_check(file):
    try:
        return make_manifold(file)
    except Exception as e:
        logging.error(f"Error in {file}: {str(e)}")
        return False
    
if __name__ == "__main__":
    os.makedirs(DATASET_SAVE_DIR.replace("test/", ""), exist_ok=True)
    os.makedirs(DATASET_SAVE_DIR, exist_ok=True)
    os.makedirs(FAILED_SAVE_DIR, exist_ok=True)
    
    files = os.listdir(RAW_DATA_DIR)
    
    def init_worker():
        os.setpgrp()
    
    with ProcessPoolExecutor(initializer=init_worker) as executor:
        futures = [executor.submit(manifold_check, file) for file in files]
        for future in tqdm(as_completed(futures), total=len(futures)):
            try:
                future.result()
            except Exception as e:
                logging.error(f"Future execution failed: {str(e)}")