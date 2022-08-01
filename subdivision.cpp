#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

#define _USE_MATH_DEFINES

// Defining a 3D vector as 3 float
struct vector3 {
    float x,y,z;
};

// Function to find the next edge after a given edge
int next(int edge);

// Returns true if two vectors are the same
bool compareTwoVectors(vector3 first, vector3 second);

// computes an edge vertex given an edge
vector3 computeEdgeVertex(int e);

// computes the new position of old vertices
vector3 computeNewPosition(vector3 v, int index);

// Will hold data from the directed edge file
std::vector<vector3> vertices;
std::vector<unsigned int> firstDirectedEdge;
std::vector<unsigned int> faceVertices;
std::vector<unsigned int> otherHalf;

int main(int argc, char *argv[])
{
    // Exit program if no model is given, or if more than 1 model is given
    if(argc != 2){
        std::cout << "Usage: " << argv[0] << "  PATH_TO/<modelname>.diredge\n";
        return 0;
    };

    // will hold the name of the given model
    std::string modelName;

    // Open given .diredge file
    std::ifstream model(argv[1]);

    if(!model.is_open()){
        std::cout << "Could not open file\n";
    }else{
        // First get the model name from the parameter
        std::string stringParameter(argv[1]);
        int firstIndex = stringParameter.find_last_of("/") + 1;
        int lastIndex = stringParameter.find_last_of(".");
        int nameSize = lastIndex - firstIndex;
        modelName = stringParameter.substr(firstIndex, nameSize);

        // Store data from the file
        std::string data;

        // Go through data and add the info to the relevant arrays
        while(model.good()){

            // Take the first word of the line
            model >> data;

            if(model.eof())
                break;

            // skip comments
            if(data == "#")
                std::getline(model,data);

            // holds skippable values
            int skip;

            // read in vertices
            if(data == "Vertex"){
                vector3 vec;
                model >> skip >> vec.x >> vec.y >> vec.z;
                vertices.push_back(vec);
            }

            // read in first directed edges
            if(data == "FirstDirectedEdge"){
                unsigned int vertexIndex, edgeIndex;
                model >> vertexIndex >> edgeIndex;
                firstDirectedEdge.push_back(vertexIndex);
                firstDirectedEdge.push_back(edgeIndex);
            }

            // read in faces
            if(data == "Face"){
                unsigned int a,b,c;
                model >> skip; // discard the index
                model >> a >> b >> c;
                faceVertices.push_back(a);
                faceVertices.push_back(b);
                faceVertices.push_back(c);
            }

            // read in other halves
            if(data == "OtherHalf"){
                unsigned int a,b;
                model >> a >> b;
                otherHalf.push_back(a);
                otherHalf.push_back(b);
            }

            // read in normals
            if(data == "Normal"){
                double a,b,c,d;
                model >> a >> b >> c >> d;
            }
        }
    }
    model.close();

    // new arrays which copy old data
    // this is to preserve old data when new data is updated
    std::vector<vector3> new_vertices = vertices;
    std::vector<unsigned int> new_faceVertices = faceVertices;

    // next available id is the size of the vertices array
    int nextEdgeId = new_vertices.size();

    // compute new sub triangles
    // go over each face of the mesh
    for(int i = 0; i < faceVertices.size(); i+=3){
        // get the vertices of this face
        vector3 v0 = vertices[faceVertices[i]];
        vector3 v1 = vertices[faceVertices[i+1]];
        vector3 v2 = vertices[faceVertices[i+2]];

        // compute new edge vertices
        vector3 v01 = computeEdgeVertex(i);
        vector3 v12 = computeEdgeVertex(i+1);
        vector3 v20 = computeEdgeVertex(i+2);

        // construct new faces
        int v01Index, v12Index, v20Index;
        v01Index = v12Index = v20Index = -1;

        // check if the edge vertices already exist
        for(int j = 0; j < new_vertices.size(); j++){
            if(compareTwoVectors(v01, new_vertices[j])){
                v01Index = j; // if this edge vertex already exists, then get set its index to the index at which it was found
            }
            if(compareTwoVectors(v12, new_vertices[j])){
                v12Index = j;
            }
            if(compareTwoVectors(v20, new_vertices[j])){
                v20Index = j;
            }
        }

        // the following if statements check to see if the edge vertices were NOT found
        // if they were not found, then give them the next edge ID, and add them to the new vertices array
        if(v01Index == -1){
            v01Index = nextEdgeId++;
            new_vertices.push_back(v01);
        }

        if(v12Index == -1){
            v12Index = nextEdgeId++;
            new_vertices.push_back(v12);
        }

        if(v20Index == -1){
            v20Index = nextEdgeId++;
            new_vertices.push_back(v20);
        }

        // now construct all of the faces by using the edge id's found, and the loop index as ID's for original vertices
        // firstly, add new faces to the new_faceVertices array
        //bottom left
        new_faceVertices.push_back(faceVertices[i]); // v0
        new_faceVertices.push_back(v01Index); // v01
        new_faceVertices.push_back(v20Index); // v20

        //bottom right
        new_faceVertices.push_back(v01Index); // v01
        new_faceVertices.push_back(faceVertices[i+1]); // v1
        new_faceVertices.push_back(v12Index); // v12

        // top
        new_faceVertices.push_back(v20Index); // v20
        new_faceVertices.push_back(v12Index); // v12
        new_faceVertices.push_back(faceVertices[i+2]); // v2

        // and secondly, update original face index to the index of the central face
        new_faceVertices[i] = v01Index;
        new_faceVertices[i+1] = v12Index;
        new_faceVertices[i+2] = v20Index;
    }

    // recalculate original vertices
    int numOldVertices = vertices.size();
    for(int i = 0; i < numOldVertices; i++){
        new_vertices[i] = computeNewPosition(new_vertices[i], i);
    }


    // now create directed edge structure
    std::vector<unsigned int> new_firstDirectedEdge;
    std::vector<unsigned int> new_otherHalf;

    for(int i = 0; i < new_faceVertices.size(); i++){
        // We know the first half of the edge but not the other half
        int firstHalf = i;
        int secondHalf = -1;

        // Get the face index values for the first half
        int firstHalfV1 = new_faceVertices[i];
        int firstHalfV2 = new_faceVertices[next(i)];

        // Go through the edges to find the other half
        for(int j = 0; j < new_faceVertices.size(); j++){
            if(new_faceVertices[j] == firstHalfV2 && new_faceVertices[next(j)] == firstHalfV1){
                    secondHalf = j;
            }
        }

        // Add first half to the other half array
        new_otherHalf.push_back(firstHalf);

        // Check if second half was found
        if(secondHalf != -1){
            // If it was found then add it to the array
            new_otherHalf.push_back(secondHalf);
        }else {
            // If it was not found then terminate the program
            std::cout << "Failed to find other half at edge: " << firstHalf << std::endl;
            return 0;
        }
    }

    // Find all the directed edges
    for(int i = 0; i < new_vertices.size(); i++){
        int directedEdge;
        // Find the first instance of an outgoing edge in the faces list
        for(int j = 0; j < new_faceVertices.size(); j++){
            if(i==new_faceVertices[j]){
                directedEdge = j;
                break;
            }
        }

        new_firstDirectedEdge.push_back(i);
        new_firstDirectedEdge.push_back(directedEdge);
    }

    // Write out to a file
    std::ofstream diredgeFile;
    std::string filePath = modelName + ".diredge";
    diredgeFile.open(filePath.c_str());

    diredgeFile << "# University of Leeds 2020-2021\n";
    diredgeFile << "# COMP 5821M Assignment 2\n";
    diredgeFile << "# Bilal Khan\n";
    diredgeFile << "# 201051660\n";
    diredgeFile << "#\n";
    diredgeFile << "# Object name: " << modelName << "\n";
    diredgeFile << "# Vertices=" << new_vertices.size() << " Faces=" << new_faceVertices.size() / 3 << "\n";
    diredgeFile << "#\n";

    // Print all vertices to the file
    for(long vertex = 0; vertex < new_vertices.size(); vertex++)
    {
        diredgeFile << "Vertex " << vertex << " " << new_vertices[vertex].x << " " << new_vertices[vertex].y << " " << new_vertices[vertex].z << "\n";
    }

    // Print all first directed edge values to the file
    for(int i = 0; i < new_firstDirectedEdge.size(); i+=2)
    {
        diredgeFile << "FirstDirectedEdge " << new_firstDirectedEdge[i] << " " << new_firstDirectedEdge[i+1] << std::endl;
    }

    // Print all faces to the file
    int index = 0;
    for (int face = 0; face < new_faceVertices.size(); face+=3)
    {
        diredgeFile << "Face " << index << " " << new_faceVertices[face] << " " << new_faceVertices[face+1] << " " << new_faceVertices[face+2] << "\n";
        index++;
    }

    // Print other halves to the file
    for (int i = 0; i < new_otherHalf.size(); i+=2)
    {
        diredgeFile << "OtherHalf " << new_otherHalf[i] << " " << new_otherHalf[i+1] << std::endl;
    }

    diredgeFile.close();

    std::cout << "\n//////////////////////////////\n";
    std::cout << "Subdivided " << modelName << std::endl;;
    std::cout << "//////////////////////////////\n\n";

    return 0;

}

// Returns true if the two vectors are the same
bool compareTwoVectors(vector3 first, vector3 second){
    if(first.x == second.x && first.y == second.y && first.z == second.z){
        return true;
    }else{
        return false;
    }
}

vector3 computeEdgeVertex(int e){
    // Get the two vertices on this edge
    vector3 v1 = vertices[faceVertices[e]];
    vector3 v2 = vertices[faceVertices[next(e)]];

    // get the other two vertices
    vector3 v3 = vertices[faceVertices[next(next(e))]];
    vector3 v4 = vertices[faceVertices[next(next(otherHalf[(e*2)+1]))]];

    // calculate edge vertex
    vector3 edgeVertex;
    edgeVertex.x = (3.0/8.0)*(v1.x + v2.x) + (1.0/8.0)*(v3.x + v4.x);
    edgeVertex.y = (3.0/8.0)*(v1.y + v2.y) + (1.0/8.0)*(v3.y + v4.y);
    edgeVertex.z = (3.0/8.0)*(v1.z + v2.z) + (1.0/8.0)*(v3.z + v4.z);

    return edgeVertex;
}

vector3 computeNewPosition(vector3 v, int index){
    // get all of outgoing edges from this vertex
    std::vector<unsigned int> outgoingEdges;

    // get the FDE of this vertex
    int FDE = firstDirectedEdge[(2*index) + 1];
    outgoingEdges.push_back(FDE);
    // get its other half
    int OH = otherHalf[(2 * FDE) + 1];
    // get the next edge
    int nxt = next(OH);
    outgoingEdges.push_back(nxt);

    // Get all the outgoing edges of the given vertex
    while(true){
        OH = otherHalf[(nxt*2)+1];
        nxt = next(OH);

        if(nxt == FDE){
            break;
        }else{
            outgoingEdges.push_back(nxt);
        }
    }


    // n = num neighbours
    int n = outgoingEdges.size();

    // calculate alpha
    float alpha;
    // values to make writing the formula for alpha easier
    float a = 1.0/n;
    float b = 5.0/8.0;
    float c = 3.0/8.0;
    float d = (1.0/4.0)*cos((2*M_PI) / n);

    alpha = a*(b-pow(c+d, 2));

    if(n == 3)
        alpha = 3.0/16.0;

    // calculate sum of the neighbour vertices
    vector3 sum;
    sum.x = sum.y = sum.z = 0;

    // reusable vector3 and edge integer
    vector3 neighbour;
    int outgoingEdge;

    // go through each outgoing edge of the vertex
    for(int n = 0; n < outgoingEdges.size(); n++){
        // get the neighbour vertex
        outgoingEdge = faceVertices[next(outgoingEdges[n])];
        neighbour = vertices[outgoingEdge];

        // increment the x,y,z of the sum
        sum.x = sum.x + neighbour.x;
        sum.y = sum.y + neighbour.y;
        sum.z = sum.z + neighbour.z;

    }

    vector3 newV;
    // compute the new position
    newV.x = ((1.0f - (n*alpha)) * v.x) + (alpha * sum.x);
    newV.y = ((1.0f - (n*alpha)) * v.y) + (alpha * sum.y);
    newV.z = ((1.0f - (n*alpha)) * v.z) + (alpha * sum.z);

    return newV;
}

int next(int edge){
    int face = edge / 3; // which face its on
    int which = edge % 3; // which edge of the face
    int nextWhich = (which+1) % 3; // which edge is the next one on this face

    int edgeID = face * 3 + nextWhich; // id of the next edge

    return edgeID;
}
