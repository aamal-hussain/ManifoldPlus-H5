//
// Created by Aamal Hussain on 17/12/2024.
//

#include "PXVTKDataset.h"
#include <vtkPolyDataNormals.h>
#include <vtkFeatureEdges.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkTriangle.h>
#include <vector>

PXVTKDataset::PXVTKDataset(const MatrixD &verts, const MatrixI &faces) : verts_(verts), faces_(faces) {
    const vtkNew<vtkPoints> points;
    const vtkNew<vtkCellArray> cells;

    for (int i = 0; i < verts.rows(); i++) {
        points->InsertNextPoint(
            verts(i, 0),
            verts(i, 1),
            verts(i, 2));
    }

    for (int i = 0; i < faces.rows(); i++) {
        std::vector<vtkIdType> cell_ids(faces.cols());
        for (int j = 0; j < faces.cols(); j++) {
            cell_ids[j] = static_cast<vtkIdType>(faces(i, j));
        }
        cells->InsertNextCell(faces.cols(), cell_ids.data());
    }

    mesh_->SetPoints(points);
    mesh_->SetPolys(cells);

    areas_.resize(faces.rows());
    has_cell_areas_ = false;
}

void PXVTKDataset::ComputeNormals(const bool compute_cell_normals, const bool compute_point_normals) const {
    const vtkNew<vtkPolyDataNormals> normal_generator;

    normal_generator->SetInputData(mesh_);
    normal_generator->SplittingOff();
    if (compute_point_normals) {
        normal_generator->ComputePointNormalsOn();
    }
    if (compute_cell_normals) {
        normal_generator->ComputeCellNormalsOn();
    }

    normal_generator->Update();

    mesh_->ShallowCopy(normal_generator->GetOutput());
}

void PXVTKDataset::ComputeCellAreas() {
    for (int i=0; i < mesh_->GetNumberOfCells(); i++) {
        vtkCell* vcell = mesh_->GetCell(i);
        const auto triangle = dynamic_cast<vtkTriangle*>(vcell);
        double p0[3], p1[3], p2[3];
        triangle->GetPoints()->GetPoint(0, p0);
        triangle->GetPoints()->GetPoint(1, p1);
        triangle->GetPoints()->GetPoint(2, p2);

        areas_(i) = vtkTriangle::TriangleArea(p0, p1, p2);
    }

    has_cell_areas_ = true;
}

MatrixD PXVTKDataset::GetVerts() const {
    const vtkIdType num_points = mesh_->GetNumberOfPoints();
    MatrixD verts_matrix(num_points, 3);

    for (vtkIdType i = 0; i < num_points; ++i) {
        double point[3];
        mesh_->GetPoint(i, point);
        verts_matrix(i, 0) = point[0];
        verts_matrix(i, 1) = point[1];
        verts_matrix(i, 2) = point[2];
    }

    return verts_matrix;
}


MatrixI PXVTKDataset::GetFaces() const {
    const vtkIdType num_faces = mesh_->GetNumberOfCells();
    MatrixI faces_matrix(num_faces, 3);

    for (vtkIdType i = 0; i < num_faces; i++) {
        vtkCell* cell = mesh_->GetCell(i);
        for (int j = 0; j < 3; j++) {
            faces_matrix(i, j) = static_cast<int>(cell->GetPointId(j));
        }
    }

    return faces_matrix;
}

MatrixD PXVTKDataset::GetCellNormals() const {
    vtkDataArray* normals = mesh_->GetCellData()->GetNormals();
    if (!normals)
        throw std::runtime_error("Cell normals not computed");

    const vtkIdType num_faces = mesh_->GetNumberOfCells();
    const auto normals_ptr = static_cast<const double*>(normals->GetVoidPointer(0));

    return Eigen::Map<const Eigen::MatrixXd>(normals_ptr, num_faces, 3);
}

MatrixD PXVTKDataset::GetPointNormals() const {
    vtkDataArray* normals = mesh_->GetPointData()->GetNormals();
    if (!normals) {
        throw std::runtime_error("Point normals not computed");
    }
    const vtkIdType num_points = mesh_->GetNumberOfPoints();
    const auto normals_ptr = static_cast<const double*>(normals->GetVoidPointer(0));

    return Eigen::Map<const Eigen::MatrixXd>(normals_ptr, num_points, 3);
}

Eigen::VectorXd PXVTKDataset::GetCellAreas() const {
    if (!has_cell_areas_) {
        throw std::runtime_error("Cell areas not computed");
    }

    return areas_;
}

ManifoldStatus PXVTKDataset::checkManifold() const {
    ManifoldStatus status{};

    const vtkNew<vtkFeatureEdges> feature_edges;
    feature_edges->SetInputData(mesh_);
    feature_edges->BoundaryEdgesOn();
    feature_edges->NonManifoldEdgesOn();
    feature_edges->ManifoldEdgesOff();
    feature_edges->FeatureEdgesOff();
    feature_edges->Update();

    auto* output = feature_edges->GetOutput();
    status.num_boundary_edges = output->GetNumberOfLines();
    status.num_non_manifold_edges = feature_edges->GetNonManifoldEdges();

    status.num_non_triangular_faces = 0;
    for (vtkIdType i = 0; i < mesh_->GetNumberOfCells(); i++) {
        if (mesh_->GetCell(i)->GetNumberOfPoints() != 3) {
            status.num_non_triangular_faces++;
        }
    }

    status.is_manifold = (status.num_boundary_edges == 0 &&
                         status.num_non_manifold_edges == 0 &&
                         status.num_non_triangular_faces == 0);

    return status;
}
