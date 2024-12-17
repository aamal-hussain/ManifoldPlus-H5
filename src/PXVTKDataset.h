//
// Created by Aamal Hussain on 17/12/2024.
//

#ifndef PXVTKDATASET_H
#define PXVTKDATASET_H
#include <vtkNew.h>
#include <vtkPolyData.h>

#include "types.h"


class PXVTKDataset {
public:
    PXVTKDataset(const MatrixD &verts, const MatrixI &faces);
    void ComputeNormals(bool compute_cell_normals = true, bool compute_point_normals = true);
    void ComputeCellAreas();

    bool HasCellAreas() const { return has_cell_areas_; }

    MatrixD GetVerts() const;
    MatrixI GetFaces() const;
    MatrixD GetCellNormals() const;
    MatrixD GetPointNormals() const;
    Eigen::VectorXd GetCellAreas() const;

private:
    vtkNew<vtkPolyData> mesh_;
    Eigen::VectorXd areas_;

    bool has_cell_areas_;
};

#endif //PXVTKDATASET_H
