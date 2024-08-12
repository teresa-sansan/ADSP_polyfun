#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <stdint.h>
#include <Eigen/Dense>

// int main() {
//     try {
//         std::string fileAnnot = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/4SBayesRC/new_all.tsv";
//         int numSnps = 10;
//         int numAnno = 200;

//         // Resize annoMat to hold numSnps * numAnno floats
//         std::vector<float> annoMat(numSnps * numAnno);

//         // Open the file in binary mode
//         FILE *hAnnot = fopen(fileAnnot.c_str(), "rb");

//         // Calculate the number of elements
//         int64_t n_ele = static_cast<int64_t>(numSnps) * numAnno;

//         // Read the data from the file
//         if (fread(annoMat.data(), sizeof(float), n_ele, hAnnot) != n_ele) {
//             std::cerr << "Error: read annotation file error, size is not correct: " << fileAnnot << std::endl;
//             throw std::runtime_error("Error reading file");
//         }

//         // Close the file
//         fclose(hAnnot);

//         // Output the read data (optional, for verification)
//         for (int i = 0; i < numSnps; ++i) {
//             for (int j = 0; j < numAnno; ++j) {
//                 std::cout << annoMat[i * numAnno + j] << " ";
//             }
//             std::cout << std::endl;
//         }
        
//         std::cout << "File read successfully!" << std::endl;

//     } catch (const std::exception &e) {
//         std::cerr << "Exception: " << e.what() << std::endl;
//     }

//     return 0;
// }


using namespace std;
using namespace Eigen;

class AnnoProb {
public:
    AnnoProb(string fileAnnot, int numAnno, const VectorXf &Pi, MatrixXf &snpPi, bool bOutDetail, double initSS);

private:
    bool bOutDetail;
    int numDist;
    int numComp;
    int numAnno;
    int numSnps;
    MatrixXf annoMat;
    VectorXf annoSD;
    VectorXf vg_enrich_qt;
    VectorXf bAnnoBinary;
    VectorXf invFrac;
    MatrixXf snpP;
    MatrixXf alpha;
    VectorXf sigmaSq;
    VectorXf ssq;
};

AnnoProb::AnnoProb(string fileAnnot, int numAnno, const VectorXf &Pi, MatrixXf &snpPi, bool bOutDetail, double initSS) {
    this->bOutDetail = bOutDetail;

    numDist = Pi.size();
    numComp = Pi.size() - 1;
    this->numAnno = numAnno;
    numSnps = snpPi.rows();

    // read annot
    annoMat.resize(numSnps, numAnno);
    FILE *hAnnot = fopen(fileAnnot.c_str(), "rb");
    if(!hAnnot){
        cerr << "Error: can't open the annotation file: " << fileAnnot << endl;
        throw runtime_error("Error");
    }
    int64_t n_ele = (int64_t) numSnps * numAnno;
    if(fread(annoMat.data(), sizeof(float), n_ele, hAnnot) != n_ele){
        cerr << "Error: read annotation file error, size is not correct: " << fileAnnot << endl;
        throw runtime_error("Error");
    }
    fclose(hAnnot);

    annoSD.resize(numAnno);
    vg_enrich_qt.resize(numAnno); 
    bAnnoBinary.resize(numAnno);

    int numBinary = 0;
    int numQT = 0;
    for(int i = 0; i < numAnno; i++){
        bool bBinary = true;
        for(int j = 0; j < numSnps; j++){
            float temp = annoMat(j, i);
            if(abs(temp) > 1e-6 && abs(temp - 1.0) > 1e-6){
                bBinary = false;
                break;
            }
        }
        if(bBinary){
            annoSD[i] = 1.0;
            numBinary++;
        } else {
            annoSD[i] = sqrt((annoMat.col(i).array() - annoMat.col(i).mean()).square().sum() / numSnps);
            numQT++;
        }
        bAnnoBinary[i] = bBinary;
    }

    cout << "Total annotation category (including intercept): " << numAnno << "." << endl;
    cout << "Number of binary annotation: " << numBinary << ", quantitative annotation: " << numQT << "." << endl;

    invFrac = numSnps * 1.0 / annoMat.colwise().sum().array();

    if (numSnps != snpPi.rows()) {
        throw runtime_error("annoMat should have same rows with snpPi");
    }
    if (snpPi.cols() != Pi.size()) {
        throw runtime_error("snpPi should have same columns with Pi");
    }

    // init the snpP
    snpP.resize(numSnps, numComp);
    // init alpha and snpPi
    alpha = MatrixXf::Zero(numAnno, numComp);
    
    // Init sigmaSS
    sigmaSq = VectorXf::Ones(numComp) * (float)initSS;
    ssq.resize(numComp);
}

int main() {
    // Example usage of the AnnoProb constructor
    string fileAnnot = "path/to/annotation/file";
    int numAnno = 10;
    VectorXf Pi(5);
    Pi << 0.2, 0.2, 0.2, 0.2, 0.2;
    MatrixXf snpPi(100, 5);
    snpPi.setConstant(0.2);
    bool bOutDetail = true;
    double initSS = 1.0;

    try {
        AnnoProb annoProb(fileAnnot, numAnno, Pi, snpPi, bOutDetail, initSS);
    } catch (const runtime_error &e) {
        cerr << e.what() << endl;
        return 1;
    }

    return 0;
}