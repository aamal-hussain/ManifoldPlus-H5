//
// Created by Aamal Hussain on 17/12/2024.
//

#ifndef PXVTKDATASET_H
#define PXVTKDATASET_H
#include <vtkNew.h>
#include <vtkPolyData.h>

#include "types.h"

struct ManifoldStatus
{
    bool is_manifold{false};
    vtkIdType num_boundary_edges{0};
    vtkIdType num_non_manifold_edges{0};
    vtkIdType num_non_triangular_faces{0};
};

class PXVTKDataset
{
public:
    PXVTKDataset(const MatrixD &verts, const MatrixI &faces);
    void ComputeNormals(bool compute_cell_normals = true, bool compute_point_normals = true) const;
    void ComputeCellAreas();

    [[nodiscard]] bool HasCellAreas() const { return has_cell_areas_; }

    MatrixD GetVerts() const;
    MatrixI GetFaces() const;
    MatrixD GetCellNormals() const;
    MatrixD GetPointNormals() const;
    Eigen::VectorXd GetCellAreas() const;

    ManifoldStatus checkManifold() const;

private:
    vtkSmartPointer<vtkPolyData> mesh_{vtkSmartPointer<vtkPolyData>::New()};
    Eigen::VectorXd areas_;

    MatrixD verts_;
    MatrixI faces_;

    bool has_cell_areas_;
};

#endif // PXVTKDATASET_H
