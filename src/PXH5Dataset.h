//
// Created by Aamal Hussain on 16/12/2024.
//

#ifndef PXH5DATASET_H
#define PXH5DATASET_H

#include <string>

#include <highfive/H5File.hpp>
#include "types.h"
#include <filesystem>
#include "PXVTKDataset.h"

using HighFive::File;
using HighFive::Group;
using HighFive::ObjectType;
namespace fs = std::filesystem;

struct DatasetParameters
{
    std::string verts;
    std::string faces;
    std::string verts_normals;
    std::string faces_normals;
    std::string areas;

    DatasetParameters(
        std::string verts,
        std::string faces,
        std::string verts_normals,
        std::string faces_normals,
        std::string areas) : verts(std::move(verts)),
                             faces(std::move(faces)),
                             verts_normals(std::move(verts_normals)),
                             faces_normals(std::move(faces_normals)),
                             areas(std::move(areas)) {}
};

class PXH5Dataset
{

public:
    explicit PXH5Dataset(DatasetParameters dataset_parameters) : dataset_parameters_(std::move(dataset_parameters)) {}
    PXH5Dataset(
        std::string verts,
        std::string faces,
        std::string verts_normals,
        std::string faces_normals,
        std::string areas) : dataset_parameters_(DatasetParameters(std::move(verts),
                                                                   std::move(faces),
                                                                   std::move(verts_normals),
                                                                   std::move(faces_normals),
                                                                   std::move(areas))) {}

    static void PrintGroupStructure(const Group &group, const std::string &prefix = "");

    template <typename T>
    T LoadMatrixDataset(const File &file, const std::string &path);

    void VertsAndFacesFromH5(const fs::path &filepath);
    void VertsAndFacesToH5(const fs::path &filepath, const MatrixD &out_verts, const MatrixI &out_faces) const;
    void PXVTKToH5(
        const fs::path &filepath,
        const PXVTKDataset &polydata,
        bool has_point_normals = true,
        bool has_cell_normals = true,
        bool has_areas = true) const;

    void ManifoldDataToH5(
        const fs::path &filepath,
        const MatrixD &out_verts,
        const MatrixI &out_faces,
        const MatrixD &out_verts_normals,
        const MatrixD &out_faces_normals) const;

    MatrixD GetVerts() const { return verts_; }
    MatrixI GetFaces() const { return faces_; }
    long GetNumberOfVerts() const { return verts_.rows(); }
    long GetNumberOfFaces() const { return faces_.rows(); }

private:
    DatasetParameters dataset_parameters_;
    MatrixD verts_;
    MatrixI faces_;
};

#endif // PXH5DATASET_H
