//
// Created by Aamal Hussain on 17/12/2024.
//

#include "PXVTKDataset.h"
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkTriangle.h>

PXVTKDataset::PXVTKDataset(const MatrixD &verts, const MatrixI &faces) {
    const vtkNew<vtkPoints> points;
    const vtkNew<vtkCellArray> cells;

    for (int i = 0; i < verts.rows(); i++) {
        points->InsertNextPoint(
            verts(i, 0),
            verts(i, 1),
            verts(i, 2));
    }

    for (int i = 0; i < faces.rows(); i++) {
        cells->InsertNextCell(3, faces.row(i).data());
    }

    mesh_->SetPoints(points);
    mesh_->SetPolys(cells);

    areas_.resize(faces.rows());
    has_cell_areas_ = false;
}

void PXVTKDataset::ComputeNormals(bool compute_cell_normals, bool compute_point_normals) {
    const vtkNew<vtkPolyDataNormals> normal_generator;

    normal_generator->SetInputData(mesh_);
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
    Eigen::MatrixXd verts(num_points, 3);
    for (vtkIdType i; i < num_points; i++) {
        double point[3];
        mesh_->GetPoint(i, point);
        verts(i, 0) = point[0];
        verts(i, 1) = point[1];
        verts(i, 2) = point[2];
    }

    return verts;
}


MatrixI PXVTKDataset::GetFaces() const {
    const vtkIdType num_faces = mesh_->GetNumberOfCells();
    Eigen::MatrixXi faces(num_faces, 3);
    for (vtkIdType i = 0; i < num_faces; i++) {
        vtkNew<vtkIdList> face;
        mesh_->GetCellPoints(i, face);
        for (int j = 0; j < 3; ++j) {
            faces(i, j) = face->GetId(j);
        }
    }

    return faces;
}

MatrixD PXVTKDataset::GetCellNormals() const {
    const vtkDataArray* normals {mesh_->GetCellData()->GetNormals()};
    if (!normals)
        throw std::runtime_error("Cell normals not computed");

    const vtkIdType num_faces = mesh_->GetNumberOfCells();
    Eigen::MatrixXd normal_matrix(num_faces, 3);

    for (vtkIdType i = 0; i < num_faces; i++) {
        double normal[3];
        normals->GetTuple(i, normal);

        normal_matrix(i, 0) = normal[0];
        normal_matrix(i, 1) = normal[1];
        normal_matrix(i, 2) = normal[2];
    }
    return normal_matrix;
}

MatrixD PXVTKDataset::GetPointNormals() const {
    const vtkDataArray* normals {mesh_->GetPointData()->GetNormals()};
    if (!normals) {
        throw std::runtime_error("Point normals not computed");
    }

    const vtkIdType num_points = mesh_->GetNumberOfPoints();
    Eigen::MatrixXd normal_matrix(num_points, 3);

    for (vtkIdType i = 0; i < num_points; i++) {
        double normal[3];
        normals->GetTuple(i, normal);

        normal_matrix(i, 0) = normal[0];
        normal_matrix(i, 1) = normal[1];
        normal_matrix(i, 2) = normal[2];
    }

    return normal_matrix;
}

Eigen::VectorXd PXVTKDataset::GetCellAreas() const {
    if (!has_cell_areas_) {
        throw std::runtime_error("Cell areas not computed");
    }

    return areas_;
}
