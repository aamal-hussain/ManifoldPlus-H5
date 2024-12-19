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

PXVTKDataset::PXVTKDataset(const MatrixD &verts, const MatrixI &faces) : verts_(verts), faces_(faces)
{
    if (verts.rows() == 0 || verts.cols() != 3)
        throw std::runtime_error("Invalid vertices matrix");

    if (faces.rows() == 0 || faces.cols() != 3)
        throw std::runtime_error("Invalid faces matrix");

    const vtkNew<vtkPoints> points;
    const vtkNew<vtkCellArray> cells;

    points->SetNumberOfPoints(verts.rows());
    for (Eigen::Index i = 0; i < verts.rows(); i++)
        points->SetPoint(
            i,
            verts(i, 0),
            verts(i, 1),
            verts(i, 2));

    cells->AllocateEstimate(faces.rows(), 3);
    for (Eigen::Index i = 0; i < faces.rows(); i++)
    {
        const vtkIdType cell_points[3] = {
            static_cast<vtkIdType>(faces(i, 0)),
            static_cast<vtkIdType>(faces(i, 1)),
            static_cast<vtkIdType>(faces(i, 2))};
        cells->InsertNextCell(3, cell_points);
    }

    mesh_->SetPoints(points);
    mesh_->SetPolys(cells);

    areas_.resize(faces.rows());
    has_cell_areas_ = false;
}

void PXVTKDataset::ComputeNormals(const bool compute_cell_normals, const bool compute_point_normals) const
{
    const vtkNew<vtkPolyDataNormals> normal_generator;

    if (!mesh_)
        throw std::runtime_error("Mesh is not initialized");

    normal_generator->SetInputData(mesh_);
    normal_generator->SplittingOff();
    if (compute_point_normals)
    {
        normal_generator->ComputePointNormalsOn();
    }
    if (compute_cell_normals)
    {
        normal_generator->ComputeCellNormalsOn();
    }

    normal_generator->Update();

    mesh_->ShallowCopy(normal_generator->GetOutput());
}

void PXVTKDataset::ComputeCellAreas()
{
    if (!mesh_)
        throw std::runtime_error("Mesh is not initialized");

    for (int i = 0; i < mesh_->GetNumberOfCells(); i++)
    {
        vtkCell *vcell = mesh_->GetCell(i);
        if (vcell->GetCellType() != VTK_TRIANGLE)
            throw std::runtime_error("Only triangular faces are supported");

        vtkTriangle *triangle = static_cast<vtkTriangle *>(vcell);
        double p0[3], p1[3], p2[3];
        triangle->GetPoints()->GetPoint(0, p0);
        triangle->GetPoints()->GetPoint(1, p1);
        triangle->GetPoints()->GetPoint(2, p2);

        areas_(i) = vtkTriangle::TriangleArea(p0, p1, p2);
    }

    has_cell_areas_ = true;
}

MatrixD PXVTKDataset::GetVerts() const
{
    if (!mesh_)
        throw std::runtime_error("Mesh is not initialized");

    const vtkIdType num_verts{mesh_->GetNumberOfPoints()};
    MatrixD verts_matrix(num_verts, 3);

    vtkPoints *points{mesh_->GetPoints()};
    if (!points)
        throw std::runtime_error("No points in mesh");

    double p[3];
    for (vtkIdType i{0}; i < num_verts; i++)
    {
        points->GetPoint(i, p);
        verts_matrix(i, 0) = p[0];
        verts_matrix(i, 1) = p[1];
        verts_matrix(i, 2) = p[2];
    }
    return verts_matrix;
}

MatrixI PXVTKDataset::GetFaces() const
{
    if (!mesh_)
        throw std::runtime_error("Mesh is not initialized");

    const vtkIdType num_faces{mesh_->GetNumberOfCells()};
    MatrixI faces_matrix(num_faces, 3);

    vtkIdType npts;
    const vtkIdType *pts;

    vtkCellArray *cells{mesh_->GetPolys()};
    if (!cells)
        throw std::runtime_error("No faces in mesh");

    vtkIdType cell_idx{0};
    for (cells->InitTraversal(); cells->GetNextCell(npts, pts); cell_idx++)
    {
        if (npts != 3)
            throw std::runtime_error("Only triangular faces are supported");

        for (int i = 0; i < npts; i++)
            faces_matrix(cell_idx, i) = static_cast<int>(pts[i]);
    }

    return faces_matrix;
}

MatrixD PXVTKDataset::GetCellNormals() const
{
    if (!mesh_)
        throw std::runtime_error("Mesh is not initialized");

    vtkDataArray *normals = mesh_->GetCellData()->GetNormals();
    if (!normals)
        throw std::runtime_error("Cell normals not computed");

    const vtkIdType num_faces = mesh_->GetNumberOfCells();
    MatrixD normals_matrix(num_faces, 3);

    for (vtkIdType i{0}; i < num_faces; i++)
    {
        double normal[3];
        normals->GetTuple(i, normal);
        normals_matrix(i, 0) = normal[0];
        normals_matrix(i, 1) = normal[1];
        normals_matrix(i, 2) = normal[2];
    }

    return normals_matrix;
}

MatrixD PXVTKDataset::GetPointNormals() const
{
    if (!mesh_)
        throw std::runtime_error("Mesh is not initialized");

    vtkDataArray *normals = mesh_->GetPointData()->GetNormals();
    if (!normals)
        throw std::runtime_error("Point normals not computed");

    const vtkIdType num_points = mesh_->GetNumberOfPoints();
    MatrixD normals_matrix(num_points, 3);

    if (normals->GetDataType() == VTK_DOUBLE && normals->GetNumberOfComponents() == 3)
    {
        double *data = static_cast<double *>(normals->GetVoidPointer(0));
        normals_matrix = Eigen::Map<MatrixD>(data, num_points, 3);
    }
    else
    {
        double normal[3];
        for (vtkIdType i = 0; i < num_points; i++)
        {
            normals->GetTuple(i, normal);
            normals_matrix.row(i) = Eigen::Vector3d(normal);
        }
    }

    return normals_matrix;
}

Eigen::VectorXd PXVTKDataset::GetCellAreas() const
{
    if (!has_cell_areas_)
    {
        throw std::runtime_error("Cell areas not computed");
    }

    return areas_;
}

ManifoldStatus PXVTKDataset::checkManifold() const
{
    ManifoldStatus status{};

    if (!mesh_)
        throw std::runtime_error("Mesh is not initialized");

    const vtkNew<vtkFeatureEdges> feature_edges;
    feature_edges->SetInputData(mesh_);
    feature_edges->BoundaryEdgesOn();
    feature_edges->NonManifoldEdgesOn();
    feature_edges->ManifoldEdgesOff();
    feature_edges->FeatureEdgesOff();
    feature_edges->Update();

    if (!feature_edges->GetOutput())
        throw std::runtime_error("Failed to compute feature edges");

    auto *output = feature_edges->GetOutput();
    status.num_boundary_edges = output->GetNumberOfLines();
    status.num_non_manifold_edges = feature_edges->GetNonManifoldEdges();

    status.num_non_triangular_faces = 0;
    for (vtkIdType i = 0; i < mesh_->GetNumberOfCells(); i++)
    {
        if (mesh_->GetCell(i)->GetNumberOfPoints() != 3)
        {
            status.num_non_triangular_faces++;
        }
    }

    status.is_manifold = (status.num_boundary_edges == 0 &&
                          status.num_non_manifold_edges == 0 &&
                          status.num_non_triangular_faces == 0);

    return status;
}
