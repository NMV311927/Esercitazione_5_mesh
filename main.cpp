#include <iostream>
#include <list>
#include <map>
#include "PolygonalMesh.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;

int main()
{
    PolygonalMesh mesh;

    if(!ImportMesh(mesh))
    {
        cerr << "file not found" << endl;
        return 1;
    }

    Gedim::UCDUtilities utilities;
    {
        vector<Gedim::UCDProperty<double>> cell0Ds_properties(1);

        cell0Ds_properties[0].Label = "Marker";
        cell0Ds_properties[0].UnitLabel = "-";
        cell0Ds_properties[0].NumComponents = 1;

        vector<double> cell0Ds_marker(mesh.NumCell0Ds, 0.0);
        for(const auto &m : mesh.MarkerCell0Ds)
            for(const unsigned int id: m.second)
                cell0Ds_marker.at(id) = m.first;

        cell0Ds_properties[0].Data = cell0Ds_marker.data();

        utilities.ExportPoints("./Cell0Ds.inp",
                               mesh.Cell0DsCoordinates,
                               cell0Ds_properties);
    }

    {

        vector<Gedim::UCDProperty<double>> cell1Ds_properties(1);

        cell1Ds_properties[0].Label = "Marker";
        cell1Ds_properties[0].UnitLabel = "-";
        cell1Ds_properties[0].NumComponents = 1;

        vector<double> cell1Ds_marker(mesh.NumCell1Ds, 0.0);
        for(const auto &m : mesh.MarkerCell1Ds)
            for(const unsigned int id: m.second)
                cell1Ds_marker.at(id) = m.first;

        cell1Ds_properties[0].Data = cell1Ds_marker.data();

        utilities.ExportSegments("./Cell1Ds.inp",
                                 mesh.Cell0DsCoordinates,
                                 mesh.Cell1DsExtrema,
                                 {},
                                 cell1Ds_properties);
    }


    // Markers check

    cout << endl << "Markers Check" << endl << "(If marker not printed then zero)" << endl;
    cout << "MarkerCell0Ds: " << endl;
    for(const auto& pair: mesh.MarkerCell0Ds){
        cout << pair.first << ": ";
        for(const auto& value : pair.second){
            cout << value << " ";
        }
        cout << endl;
    }
    cout << endl;
    
    cout << "MarkerCell1Ds: " << endl;
    for(const auto& pair: mesh.MarkerCell1Ds){
        cout << pair.first << ": ";
        for(const auto& value : pair.second){
            cout << value << " "; 
        }
        cout << endl;

    }
    cout << endl;


    cout << "MarkerCell2Ds: " << endl;
    for(const auto& pair: mesh.MarkerCell2Ds){
        cout << pair.first << ": ";
        for(const auto& value : pair.second){
            cout << value << " ";
        }
        cout << endl;
    }
    cout << "(Correct if none printed as all should be zero)" << endl;


    // Edge length check

    const double err = 1e-10;

    for (unsigned int i = 0; i < mesh.NumCell1Ds; i++){
       
        unsigned int id0 = mesh.Cell1DsExtrema(0, i);
        unsigned int id1 = mesh.Cell1DsExtrema(1, i);
        Vector3d v0 = mesh.Cell0DsCoordinates.col(id0);
        Vector3d v1 = mesh.Cell0DsCoordinates.col(id1);

        double length = (v1 - v0).norm();

        if (length < err){
            cout << "Careful! Edge " << i << " has zero length" << endl;
        }
            
        }
    


    // Polygon area check

	for (unsigned int i = 0; i < mesh.NumCell2Ds; i++){
        vector<unsigned int> vertices = mesh.Cell2DsVertices[i];

        if (vertices.size() < 3){
            double area = 0.0;
            for (unsigned int j = 0; j < vertices.size(); ++j){
                unsigned int k = (j + 1) % vertices.size();
                Vector2d v0 = mesh.Cell0DsCoordinates.col(vertices[j]);
                Vector2d v1 = mesh.Cell0DsCoordinates.col(vertices[k]);

                area += v0(0) * v1(1) - v0(1) * v1(0);
                area = abs(area) / 2.0;
            }

            if (area < err){
                cout << "Careful! Polygon " << i << " has zero area" << endl;
        
            }
        }
    }


    return 0;
}
