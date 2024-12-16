#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include "Manifold.h"
#include "Parser.h"
#include "types.h"
#include "PXH5Dataset.h"

#include <filesystem>

namespace fs = std::__fs::filesystem;

void processFile(const fs::path &source, const fs::path &target, const int depth, PXH5Dataset &dataset)
{

	dataset.ReadH5(source);
	Manifold manifold;
	std::cout << "Number of verts: " << dataset.GetNumberOfVerts() << std::endl;
	std::cout << "Number of faces: " << dataset.GetNumberOfFaces() << std::endl;
	MatrixD out_verts;
	MatrixI out_faces;

	manifold.ProcessManifold(dataset.GetInVerts(), dataset.GetInFaces(), depth, &out_verts, &out_faces);

	dataset.WriteH5(target, out_verts, out_faces);
}

int main(const int argc, char **argv)
{
	Parser parser;
	parser.AddArgument("source", "../examples/input.h5");
	parser.AddArgument("target", "../examples/output.h5");
	parser.AddArgument("depth", "8");

	parser.AddArgument("verts", "mesh.verts");
	parser.AddArgument("faces", "mesh.faces");
	parser.AddArgument("verts_normals", "mesh.verts_normals");
	parser.AddArgument("faces_normals", "mesh.faces_normals");
	parser.AddArgument("areas", "mesh.areas");

	parser.ParseArgument(argc, argv);
	parser.Log();

	const DatasetParameters dataset_parameters(
		parser["verts"],
		parser["faces"],
		parser["verts_normals"],
		parser["faces_normals"],
		parser["areas"]);

	PXH5Dataset dataset(dataset_parameters);

	const int depth = std::stoi(parser["depth"]);

	fs::path source_filepath = parser["source"];
	fs::path target_filepath = parser["target"];

	std::cout << "Processing file " << source_filepath << std::endl;
	try
	{
		processFile(source_filepath, target_filepath, depth, dataset);
		std::cout << "Processed: " << source_filepath << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << "Failed to process " << source_filepath << ": " << e.what() << std::endl;
	}

	return 0;
}