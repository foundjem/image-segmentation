// Image Processing
// Dr. Layachi Bentabet
// Assignment 3

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include <algorithm>

using namespace std;

class Image {
  private:
    //instance vars
    string name;
    string type;
    string comments;
    int width;
    int height;
    int maxPossibleColourLevel; //can be colour or greyscale
    int maxColourLevel;   //can be colour or grey level
    int minColourLevel;   //can be colour or grey level
    vector<int> pixels;
    vector<int> original;
    vector<int> last;
    vector<int> allPositions;
    vector< vector<int> > computedFeatures;
    vector< vector<int> > objects;
    vector<int> redChannel;
    vector<int> greenChannel;
    vector<int> blueChannel;
    
  public:
    //constructor
    Image (char* filename) {
        name = filename;
        ifstream imageFile(filename);
        string line;
        pixels = vector<int>();
        original = vector<int>();
        last = vector<int>();
        allPositions = vector<int>();
        objects = vector< vector<int> >();
        if (imageFile.is_open()) {
            cout << " |    Successfully opened " << name << endl;
            getline(imageFile, type);
            if (type=="P2" || type=="P3") {
                cout << " |      Type: " << (type=="P2"?"Greyscale":"Colour") << endl;
                while (imageFile.peek() == '#') {
                    getline(imageFile, line);
                    comments = comments + line;
                }
                cout << " |      Comments: " << (comments==""?"(none)":comments) << endl;
                imageFile >> width;
                cout << " |      Width: " << width << endl;
                imageFile >> height;
                cout << " |      Height: " << height << endl;
                imageFile >> maxPossibleColourLevel;
                cout << " |      Maximum possible "<< (type=="P2"?"grey ":"colour ") << "level: " << maxPossibleColourLevel << endl;
                maxColourLevel = 0;
                minColourLevel = maxPossibleColourLevel;
                int temp;
                if (type=="P2") {
                    for (int i = 0; i<width*height; ++i) {
                        imageFile >> temp;
                        pixels.push_back(temp);
                        original.push_back(temp);
                        last.push_back(temp);
                        allPositions.push_back(i);

                        if (pixels[i]>maxColourLevel)
                            maxColourLevel = pixels[i];
                        if (pixels[i]<minColourLevel)
                            minColourLevel = pixels[i];
                    }
                    cout << " |      Actual maximum grey level: " << maxColourLevel << endl;
                    cout << " |      Actual minimum grey level: " << minColourLevel << endl;
                }
                else {
                    for (int i = 0; i<width*height; ++i) {
                        int sum;
                        imageFile >> temp;
                        redChannel.push_back(temp);
                        sum = temp-temp%32;
                        imageFile >> temp;
                        greenChannel.push_back(temp);
                        sum += (temp-temp%32)/8;
                        imageFile >> temp;
                        blueChannel.push_back(temp);
                        sum += (temp-temp%64)/32;
                        pixels.push_back(sum);
                        original.push_back(sum);
                        last.push_back(sum);
                        allPositions.push_back(i);

                        if (pixels[i]>maxColourLevel)
                            maxColourLevel = pixels[i];
                        if (pixels[i]<minColourLevel)
                            minColourLevel = pixels[i];
                    }
                }
            }
            else cerr << " |     Sorry, only P1 or P2 can be read.           |" << endl;
            // IS THERE A WAY TO CALL IMAGEMAGICK FROM HERE?
            imageFile.close();
        }
        else cerr << " |      Unable to open file " << endl;  
    }
    

    //======================================================================================
    //get methods

    int getPixel (int position) {
      return pixels[position];
    }
    int getX (int position) {
        return position%width;
    }
    int getY (int position) {
        return (position-(position%width))/width;
    }
    vector<int> getAllPixels () {
        return allPositions;
    }
    vector<int> getNeighbours (int position, int neighbourhoodSize) {
        vector<int> neighboursPositions = vector<int>();
        int Xposition = getX(position);
        int Yposition = getY(position);
        for (int i=0; i<((2*neighbourhoodSize)+1); ++i) {
            for (int j=0; j<((2*neighbourhoodSize)+1); j++) {
                if ( (Xposition-neighbourhoodSize+i>=0) && (Xposition-neighbourhoodSize+i<width) &&
                     (Yposition-neighbourhoodSize+j>=0) && (Yposition-neighbourhoodSize+j<height) ) {
                    neighboursPositions.push_back( (position-neighbourhoodSize)-(neighbourhoodSize*width)+i+(j*width) );
                }
            }
        }
        return neighboursPositions;
    }


    //======================================================================================
    //control methods

    void revert () {
        this->pixels = this->original;
        cout << " |      Image reverted to original" << endl;
    }
    void undo () {
        this->pixels = this->last;
        cout << " |      Last change undone" << endl;
    }
    void outputToFile (char* filename, vector<int> positions) {
        ofstream imageFile(filename);
        if (imageFile.is_open()) {
            imageFile << type << endl; 
            imageFile << "# created from " << name << endl; 
            imageFile << width << " " << height << endl;
            imageFile << maxPossibleColourLevel << endl;    
            for (int i=0; i<positions.size(); ++i) {
                imageFile << pixels[positions[i]] << " ";
            }
            imageFile << endl;
            imageFile.close();
            cout << " |      New file called " <<filename << " created" <<endl;
        }
        else cerr << " |      Unable to create file" << endl;
    }
 

    //======================================================================================
    //compute histograms

    vector<int> greyscaleHistogram () {
        vector<int> histogram = vector<int>(maxPossibleColourLevel);
        for (int i=0; i<pixels.size(); ++i) {
            ++histogram[pixels[i]];
        }
        computedFeatures.push_back(histogram);
        cout << " |      Histogram computed " << endl;
        return histogram;
    }
    vector<int> colourHistogram () {
        cout << "colour?" << endl;
        vector<int> histogram = greyscaleHistogram();
       return histogram;
    }
    vector<int> textureHistogram () {
        vector<int> histogram = vector<int>();
        cout << " |      Texture histogram computed " << endl;
        return histogram;
    }


    //======================================================================================
    //clustering

    void kmeans (int k, int numberOfClusters, vector<int> histogram) {
        
        cout << " |      Kmeans clustering begun" << endl;
        //initializations
        vector< vector<float> > clusterValues = vector< vector<float> >(k);
        vector<float> meansOld = initRandom(k);
        vector<float> meansNew = vector<float>(k);
        vector<int> closestMeans = vector<int>();
        vector<int> closestMeansNew = vector<int>();
        bool somePointMoves = true;
        int iterations = 0;
        int maxIterations = 60;

        //until means don't move much, assign pixels to a cluster and recalculate mean
        while (somePointMoves == true && iterations < maxIterations) {
            closestMeans = getClosestMeans(k, meansOld);
            for (int i=0; i<pixels.size(); ++i) {
                //cout << pixels[i] << endl;
                //cout << closestMeans[pixels[i]] << endl;
                clusterValues[closestMeans[pixels[i]]].push_back(pixels[i]);
            }
            cout << "here?" <<endl;
            meansNew = computeMeans(clusterValues, histogram);
            closestMeansNew = getClosestMeans(k, meansNew);
            if (different(closestMeans, closestMeansNew) == false) {
                cout << "here1?" << endl;
                somePointMoves = false;
            }
            else {
                for (int i=0; i<k; ++i) {
                    meansOld[i] = meansNew[i];
                }
            }
            ++iterations;
        }
        for (int i=0; i<pixels.size(); ++i) {
            pixels[i] = (closestMeans[pixels[i]])*maxPossibleColourLevel/k;
        }
        cout << " |      Kmeans clustering finished" << endl;
    }
    //set means to the intervals at (L/k)*i (i [0->k])
    vector<float> initIntervals (int k) {
        vector<float> means = vector<float>(k);
        for (int i=0; i<k; ++i) {
            means[i] = (maxPossibleColourLevel/k)*i;
        }
        cout << " |        Means initialized to intervals" << endl;
        return means;
    }
    //set means to random pixels
    vector<float> initRandom (int k) {
        vector<float> means = vector<float>(k);
        for (int i=0; i<k; ++i) {
            means[i] = pixels[rand()%(width*height)];
            cout << means[i] << " ";
        }
        cout << endl;
        cout << " |        Means initialized to random pixels" << endl;
        return means;
    }
    //faster assignment for larger images (L*k + n where L is usually 255)
    vector<int> getClosestMeans (int k, vector<float> means) {
        vector<int> closestMeans = vector<int>();
        float smallestDistance;
        int meanAssigned = 0;
        for (int i=0; i<maxPossibleColourLevel+1; ++i) {
            smallestDistance = maxPossibleColourLevel+5;
            for (int j=0; j<means.size(); ++j) {
                if ( sqrt(pow(i - means[j] , 2)) < smallestDistance ) {
                    smallestDistance = sqrt(pow(i - means[j] , 2));
                    meanAssigned = j;
                }
            }
            closestMeans.push_back(meanAssigned);
        }
        cout << " |        Closest means found" << endl;
        return closestMeans;
    }
    vector<float> computeMeans (vector< vector<float> > clusterValues, vector<int> histogram) {
        vector<float> means = vector<float>(clusterValues.size());
        float histogramSumForMean;
        for (int i=0; i<clusterValues.size(); ++i) { //for each mean
            histogramSumForMean = 0;
            for (int j=0; j<clusterValues[i].size(); ++j) { //for all the pixels belonging to that cluster
                histogramSumForMean += histogram[pixels[i]];
            }
            if(histogramSumForMean!=0)
            for (int j=0; j<clusterValues[i].size(); ++j) { //for all the pixels belonging to that cluster
                means[i] += clusterValues[i][j]*histogram[pixels[i]]/histogramSumForMean;
            }
            cout << means[i] << " ";
        }
        cout << " |        Means computed" << endl;
        return means;
    }
    bool different (vector<int> firstVector, vector<int> secondVector) {
        if (firstVector.size() != secondVector.size())
            return true;
        for (int i=0; i<firstVector.size(); ++i) {
            if (firstVector[i] != secondVector[i])
                return true;
        }
        return false;
    }
    void FCM (int k, int numberOfClusters) {

    }


    //======================================================================================
    //connected components analysis

    void makeBinary () {
        last = pixels;
        for (int i=0; i<pixels.size(); ++i) {
            if (pixels[i] > 150)
                pixels[i] = 255;
            else 
                pixels[i] = 0;
        }
        cout << " |      Image made binary (255=white, 0=black)" << endl;
    }
    void findComponents () {
        last = pixels;
        int objectCount = 0;
        for (int i=0; i<pixels.size(); ++i) {
            if (pixels[i] == 0) { //i.e. if not background and not yet checked
                objectCount += 1;
                pixels[i] = objectCount;
                objects.push_back(vector<int>());
                objects[pixels[i]-1].push_back(i);
                labelNeighbours(i);
            }
        }
        for (int i=0; i<objects.size(); ++i) {
            if (objects[i].size() < 20) {
                objects.erase(objects.begin()+i);
                --objectCount;
            }
        }
        cout << " |      " << objectCount << (objectCount==1?" object found":" objects found") << endl;
        colourComponents(objectCount);
    }
    void labelNeighbours (int position) {
        vector<int> currentNeighbours = getNeighbours(position, 1);
        for (int n=0; n<currentNeighbours.size(); ++n) {
            if (pixels[currentNeighbours[n]] == 0) {
                pixels[currentNeighbours[n]] = pixels[position];
                objects[pixels[position]-1].push_back(currentNeighbours[n]);
                labelNeighbours(currentNeighbours[n]);
            }
        }
    }
    void colourComponents (int objectCount) {
        last = pixels;
        int greyLevel = 200/objectCount;
        for (int i=0; i<pixels.size(); ++i) {
            if (pixels[i] != 255)
                pixels[i] = pixels[i]*greyLevel;
        }
        cout << " |      Objects coloured by greyscale" << endl;
    }
    void drawBoundingBoxes () {
        for (int i=0; i<objects.size(); ++i) {
            cout << " |      Object " << i+1 << ": " << endl;
            cout << " |        Area: " << objects[i].size() << endl;
            drawBoundingBox(objects[i], objects[i].size());
        }
    } 
    void drawBoundingBox (vector<int> positions, int area) {
        last = pixels;
        vector<int> box = getBoundingBox(positions);
        int Xmin = box[0];
        int Ymin = box[1];
        int Xmax = box[2];
        int Ymax = box[3]; 
        for (int i=Xmin; i<Xmax; ++i) {
            pixels[i+Ymin*width] = 0;
            pixels[i+Ymax*width] = 0;
        }
        for (int i=Ymin; i<Ymax; ++i) {
            pixels[Xmax+i*width] = 0;
            pixels[Xmin+i*width] = 0;
        }
        float rowMoment = area/(float)(Xmax-Xmin);
        float colMoment = area/(float)(Ymax-Ymin);
        cout << " |        Centroid: (" << rowMoment+Xmin << ", " << colMoment+Ymin <<")"<<endl;
        
        vector<int> objectSquare = getObjectSquare(Xmin, Ymin, Xmax, Ymax);
        float sum = 0.0;
        for (int j=0; j<(Ymax-Ymin); ++j) {
            int numberOfPixelsInThisRow = 0;
            for (int i=0; i<(Xmax-Xmin); ++i) {
                numberOfPixelsInThisRow+=objectSquare[i+(j*(Xmax-Xmin))]!=255?1:0; 
            }
            sum += pow((rowMoment-numberOfPixelsInThisRow),2);
        }
        float rowVariance = sum;
        sum = 0.0;
        for (int j=0; j<(Xmax-Xmin); ++j) {
            int numberOfPixelsInThisCol = 0;
            for (int i=0; i<(Ymax-Ymin); ++i) {
                numberOfPixelsInThisCol+=objectSquare[j+(i*(Xmax-Xmin))]!=255?1:0; 
            }
            sum += pow((colMoment-numberOfPixelsInThisCol),2);
        }
        float colVariance = sum;
        cout << " |        2nd order moment: (" << rowVariance << ", " << colVariance <<")"<<endl;
        cout << " |            sum: " << rowVariance+colVariance << "   ratio: " << rowVariance/colVariance << endl;
    }
    vector<int> getBoundingBox (vector<int> positions) {
        int Xmin = pixels.size();
        int Ymin = pixels.size();
        int Xmax = 0;
        int Ymax = 0
;        for (int i=0; i<positions.size(); ++i) {
            if (getX(positions[i]) > Xmax)
                Xmax = getX(positions[i]);
            if (getX(positions[i]) < Xmin)
                Xmin = getX(positions[i]);
            if (getY(positions[i]) > Ymax)
                Ymax = getY(positions[i]);
            if (getY(positions[i]) < Ymin)
                Ymin = getY(positions[i]);
        }
        cout <<" |        Bounding box ("<< Xmin <<","<< Ymin << "), ("<< Xmax <<","<< Ymax << ")" <<endl;
        vector<int> coordinates = vector<int>();
        coordinates.push_back(Xmin);
        coordinates.push_back(Ymin);
        coordinates.push_back(Xmax);
        coordinates.push_back(Ymax);
        return coordinates;
    }
    
    vector<int> getObjectSquare (int Xmin, int Ymin, int Xmax, int Ymax) {
        vector<int> positions = vector<int>();
        for(int i=Ymin; i<Ymax; ++i)
            for(int j=Xmin; j<Xmax; ++j)
                positions.push_back(pixels[i*width+j]);
        return positions;
    }
    /*
    void objectsToFiles () {
        for (int i=0; i<objects.size(); ++i) {
            stringstream name << "'object" << i << ".pgm";
            string filename = name.str();
            filename = ((char*)filename.c_str());
            outputToFile(filename, objects[i]);
        }
        cout << " |        Objects saved as separate files" << endl;
    }
    */
};


//============================================================================
//print things

void printWelcome () {
    cout << endl;
    cout << " |--------------------------------------------------|" << endl;
    cout << " |                                                  |" << endl;
    cout << " |  Welcome to command-line binary image analysis!  |" << endl;
    cout << " |  This program accepts basic (P2, P3) PGM or PBM  |" << endl;
    cout << " |  files, optionally, a set of preferences, and    |" << endl;
    cout << " |  performs image segmentation based on these      |" << endl;
    cout << " |  preferences.                                    |" << endl;
    cout << " |                                                  |" << endl;
}
void printFilePrompt () {
    cout << " |  Please type the name of your file (e.g img.pgm) |" << endl;
    cout << " |  and then hit ENTER.                             |" << endl;
    cout << " |                                                  |" << endl;
    cout << " |    File name: ";
}
void printPreferences (vector<string> features, string clusteringMethod, int numberOfClusters, int k) {
    cout << " |  " << endl;
    cout << " |    Segmentation preferences: " << endl;
    if (features.size() > 1) { 
        cout << " |      Features: "; 
        for (int i=0; i<features.size()-1; ++i) {
            cout << features[i] << ", ";
        } 
        cout << features[features.size()-1] << endl; 
    }
    else cout << " |      Feature: " << features[0] << endl; 
    cout << " |      Clustering method: " << clusteringMethod << endl;
    cout << " |      Number of clusters: " << numberOfClusters << endl;
    cout << " |      k: " << k << endl;
    cout << " |  " << endl;
    cout << " |    Changelog: " << endl;
}
void printExit() {
    cout << " |                                                  |" << endl;
    cout << " |  Exiting. Have a nice day!                       |" << endl;
    cout << " |                                               TM |" << endl;
    cout << " |--------------------------------------------------|" << endl;
    cout << endl;
}


//======================================================================================

int main (int argc, char *argv[]) {
    
    //default values
    string filename ="image.pgm";
    int numberOfClusters = 5;
    string clusteringMethod = "5means";
    int k = 5;
    vector<string> features = vector<string>();
    features.push_back("greyscaleHistogram");

    //welcome text
    printWelcome();
 
    //prompt for file if not given
    if (argc < 2) {
        printFilePrompt();
        cin >> filename;
    }
    //filename from commandline, everything else default
    else if (argc == 2) {
        filename = argv[1];
    }
    //everything specified in args 
    //(arg0, filename, #ofclusters, method, feature1, feature 2,...)
    else if (argc >=4) {
        filename = argv[1];
        stringstream ss(argv[2]);
        ss >> numberOfClusters;
        clusteringMethod = argv[3];
        features.clear();
        for (int i=4; i<argc; ++i) {
            features.push_back(argv[i]);
        }
    }
    else {
        cerr << " |    Incorrect number of parameters provided  " << endl;
        printExit();
        return 0;
    }
    //read k from kmeans (e.g. 3means, k=3)
    if (isdigit(clusteringMethod[0])) {
        k = ((int)clusteringMethod[0])-48;
    }

    //read image from file
    Image img = Image(((char*)filename.c_str()));
    printPreferences(features, clusteringMethod, numberOfClusters, k);
    
    //compute features
    vector<int> greyscaleHistogram = vector<int>();
    vector<int> colourHistogram = vector<int>();
    vector<int> textureHistogram = vector<int>();

    for (int i=0; i<features.size(); ++i) {
        if (features[i] == "greyscale" || features[i] == "Greyscale")
            greyscaleHistogram = img.greyscaleHistogram();
        else if (features[i] == "colour" || features[i] == "Colour")
            colourHistogram = img.colourHistogram();
        else if (features[i] == "texture" || features[i] == "Texture")
            textureHistogram = img.textureHistogram();
        else
            cerr << " |      No method found for feature \"" << features[i] << "\"" <<endl;
    }

    //cluster based on features
    if (clusteringMethod.substr(1) == "means") {
        img.kmeans(k, numberOfClusters, greyscaleHistogram);
        img.outputToFile("kmeans", img.getAllPixels());
    }
    if (clusteringMethod == "FCM" || clusteringMethod == "fcm")
        img.FCM(k, numberOfClusters);

    printExit();
    return 0;
}


//  ./main 1> cout 2> cerr
