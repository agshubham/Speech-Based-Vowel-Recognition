// 214101053_vowelRecognition.cpp : Defines the entry point for the console application.
// @author Shubham_214101053

#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

#define N 320  //sample count
#define p 12	// predictor order of speech
#define f 5  //frame Count

const char * recordingRollNumber = "214101053";
const char* vowelRecordingpath = "../214101053/";
//Tokhura weights
float weight[p] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

//Applying the raised sin window on Ci's


//calculate the  Ri's and printing
void calculateR(float *testDataArr, float *rArr, bool toPrint){
	
	int i = 0, j=0;
	for(i=0; i<=p; ++i) {
		for(j=0; j<=N-i-1; ++j) {
				rArr[i] += testDataArr[j] * testDataArr[j+i];
		}
	}

	if (toPrint) {
		printf("\n The Ri's values calculated for the test file are as \n ");
		for(i=0; i<=p; ++i) {
			printf("R[%d] = %f \n", i, rArr[i]);
		}
	}
}


// using levinson durbin Algo to 
// calculate the  Ai's and printing
void calculateA(float *rArr, float *aArr, bool toPrint) {

	// initialsing the array with zero
	float alpha[p+1][p+1] = {{0}};
	float e[p+1] = {0}, k[p+1] = {0}; 
	e[0] = rArr[0];

	int i=1, j=1;

	for(i=1; i<=p; i++) {
		float tempSum = 0;
		
		for(j=1; j<=i-1; j++){
			tempSum += alpha[i-1][j]*rArr[i-j];
		}

		k[i] = (rArr[i] - tempSum)/e[i-1];
		alpha[i][i] = k[i];

		for(j=1; j<=i-1; j++){
			alpha[i][j] = alpha[i-1][j] - k[i] * alpha[i-1][i-j];
		}
		e[i] = (1-k[i]*k[i])*e[i-1];
	}

	for(j=1; j<=p; j++){
		aArr[j] = alpha[p][j];
	}

	if (toPrint) {
		printf("\n The Ai's values calculated for the test file are as \n ");
		for(j=1; j<=p; j++){
				printf("A[%d] = %f \n", j, aArr[j]);
		}
	}
}


// Using cepstral Coefficients formula
//calculating Ci's and printing
void calculateC(float *cArr, float *aArr, double sigma, bool toPrint){
	
	int i = 0;
	cArr[0] = log10(sigma*sigma);

	// use Q = P for convenience 
	for(i=1; i<=p; i++) {
		float tempSum = 0;
		for(int j=1; j<=i-1; j++){
			tempSum += (j/(i*1.0))*cArr[j]*aArr[i-j];
		}
		cArr[i] = aArr[i] + tempSum;
	}
	
	if (toPrint) {
		printf("\n The Ci's values calculated for the test file are as \n ");
		for(i=1; i<=p; i++) {
			printf("c[%d] = %f \n", i, cArr[i]);
		}
	}
}

double raisedSin(int m) {
	return 1 + (p*1.0/2) * sin(3.14*m/12);
}


//Applying the hamming window on Si's
double hammingWindow(int m) {
	return (0.54 - 0.46 * cos(2*3.14*m/N-1));
}



void applyHamming(float s[][N]){
	for(int i=0; i<f; ++i){
		for(int j=0; j<N; ++j){
			s[i][j] *= hammingWindow(s[i][j]);
		}
	}
}


void printReferenceData(float (*cValueReference)[p+1]){
	int i=0,j=0;
	for (i = 0;i <f ; i++) {
		for (j=1;j<=p;j++) {
			printf("%f, ", cValueReference[i][j]);
		}
		printf("\n");
	}
}


// Verift and print the test file
void verifyTestFile(float *testDataArr){

	static float aArr[p+5];
	static float cArr[p+5];
	static float rArr[p+5];

	calculateR(testDataArr, rArr,true);
	calculateA(rArr, aArr,true);
	calculateC(cArr, aArr, rArr[0],true);

}

// to compute the average of Ci's
void findAverage(float (*output)[p+1], float (*input)[p+1]) {

	int i =0, j=0; 
	for(i=0; i<f; i++) {
		for(j=1; j<=p; j++) {
			output[i][j] += input[i][j];
			output[i][j] += input[i+5][j];
			output[i][j] += input[i+10][j];
			output[i][j] += input[i+15][j];
			output[i][j] += input[i+20][j];
			output[i][j] = output[i][j]/5;
		}
	}
}


// calculate the tokhura distance
double computeTokhuraDistance(float cTest[][p+1], float cRef[][p+1], float weight[]){
	
	double tokhuraDistance = 0;
	for(int i=0; i<f; ++i){
		double distance = 0;
		for(int j=1; j<=p; ++j){
			double diff = (cTest[i][j] - cRef[i][j]);
			distance += weight[j-1] * diff * diff;
		}
		tokhuraDistance += (distance/(float)(p));
	}
	return tokhuraDistance/(float)(f);
}

//Calculate the euclid distance
double computeEuclideanDistance(float cTest[][p+1], float cRef[][p+1]){
	
	double euclidDistance = 0;
	for(int i=0; i<f; ++i){
		double distance = 0;
		for(int j=1; j<=p; ++j){
			double diff = (cTest[i][j] - cRef[i][j]);
			distance += diff * diff;
		}
		euclidDistance += (distance/(float)(p));
	}
	return euclidDistance/(float)(f);
}

int _tmain(int argc, _TCHAR* argv[])
{
	char vowel[] = {'a','e','i','o','u'};

	// Array to store the testDataArray
	float testDataArr[320];

	int headerLineSize = 50 , headerCounter = 0, sampleLineSize = 50, count = 0 ;
	char *headerLine = (char *)malloc(headerLineSize);
	char *sampleLine = (char *)malloc(sampleLineSize);

	/**************************************************************
						PART 1 Data Verification
	****************************************************************/

	char testFileName[] = "../Test.txt";
	FILE *testFilePtr = fopen( testFileName,"r");

	if (testFilePtr != NULL) 
	{
		while (fgets(sampleLine, sampleLineSize, testFilePtr)) 
		{
			testDataArr[count++] = atof(sampleLine);
		}
		verifyTestFile(testDataArr);
	
	} else {
		float testDataArray[320] = {-53.397691, -61.915578, -46.428510, 8.550582, 11.647996, -8.485193, -3.839072, 7.001876, 105.344759, 178.133980, 129.349715, 75.144976, 35.652952, 55.786141, 121.606181, 165.744326, 154.129024, 96.826872, 41.073426, 48.816960, 77.468036, 33.329892, 34.878599, 22.488944, 11.647996, 21.714591, 21.714591, -3.064719, -44.105450, -58.818165, -39.459329, -3.839072, 25.586358, 28.683771, -13.905667, -50.300277, -47.977217, -24.746614, -14.680020, -3.839072, -2.290365, -18.551787, -8.485193, -13.905667, -23.197908, -24.746614, -12.356960, 7.776229, 3.904462, -34.038855, -85.146181, -76.628293, -64.238639, -82.823120, -67.336052, -58.818165, -87.469241, -68.884759, -78.951353, -169.550703, -248.534751, -210.591434, -131.607386, -78.177000, -71.207819, -105.279370, -75.853940, -55.720751, -24.746614, 49.591313, 79.791097, 54.237434, 147.934197, 210.656824, 165.744326, 185.877514, 168.841739, 140.190663, 194.395402, 219.949065, 192.846695, 144.836784, 77.468036, 68.950149, 79.016743, 55.011787, 72.821916, 59.657908, 8.550582, 27.135064, 30.232478, -3.064719, -16.228727, -10.033899, -5.387779, 24.812004, 54.237434, 51.140020, 1.581402, -64.238639, -70.433466, -24.746614, -9.259546, -0.741658, -17.003080, -65.012992, -67.336052, -67.336052, -65.012992, -72.756526, -95.987129, -60.366871, -23.972261, -17.003080, -20.100494, -54.172044, -90.566655, -65.787345, 13.971056, 40.299072, -14.680020, -67.336052, -79.725707, -93.664068, -91.341008, -41.008036, -42.556743, -82.823120, -116.894671, -177.294237, -224.529795, -216.786261, -166.453289, 
									-136.253506, -111.474197, -69.659112, -41.008036, -29.392735, -50.300277, -20.874847, 56.560494, 121.606181, 212.205531, 257.118029, 192.846695, 159.549498, 199.815876, 207.559410, 218.400358, 234.661780, 182.780101, 154.903378, 181.231394, 144.836784, 77.468036, 58.883554, 36.427305, 46.493900, 60.432261, 41.847779, -1.516012, -37.136269, -33.264502, 10.873643, 32.555538, 4.678815, 0.032695, -23.197908, -6.936486, -4.613426, -7.710839, -24.746614, -53.397691, -34.038855, -5.387779, -26.295321, -58.043811, -101.407602, -123.863851, -83.597474, -23.972261, -0.741658, -68.884759, -116.120317, -88.243594, -41.008036, 1.581402, 15.519763, -51.074630, -89.792301, -41.008036, 5.453169, 20.940237, -17.003080, -83.597474, -126.186912, -112.248550, -75.853940, -58.818165, -72.756526, -161.807169, -239.242510, -219.109321, -195.878719, -174.196823, -154.837988, -119.217731, -88.243594, -29.392735, 11.647996, 45.719546, 43.396486, 75.919330, 182.780101, 236.984840, 216.077298, 240.082254, 242.405314, 212.205531, 246.277081, 232.338719, 191.297988, 175.810920, 158.775145, 127.026655, 118.508767, 65.078382, 18.617177, 0.032695, 10.099289, 42.622133, 32.555538, -7.710839, -14.680020, -7.710839, 1.581402, 12.422350, -1.516012, 1.581402, -4.613426, -10.808253, -2.290365, -6.162132, -40.233683, -35.587562, -21.649201, -24.746614, -20.874847, -43.331096, -88.243594, -90.566655, -71.982173, -70.433466, -79.725707, -89.017948, -79.725707, -64.238639, -49.525924, -25.520968, -3.064719, -23.972261, -20.100494, 10.873643, -6.162132, -21.649201, 
									-35.587562, -54.946398, -63.464285, -50.300277, -16.228727, -55.720751, -150.966221, -222.206735, -254.729579, -233.822036, -160.258462, -127.735619, -154.837988, -173.422470, -170.325056, -105.279370, -1.516012, 76.693683, 112.313940, 105.344759, 161.872559, 243.179667, 247.825788, 236.984840, 247.051434, 226.143892, 233.887426, 293.512639, 252.471908, 147.159844, 87.534631, 84.437217, 98.375579, 108.442173, 71.273209, 7.001876, -20.874847, -29.392735, -6.162132, 9.324936, -13.131313, -13.905667, 11.647996, 10.099289, -10.033899, -14.680020, -29.392735, -7.710839, 17.068470, 32.555538, 7.776229, -41.782389, -45.654157, -54.172044, -51.848984 };
		verifyTestFile(testDataArray);
	}

	/****************************************************************
						PART 2  Reference Data Generation
	******************************************************************/
	int i=0, j=0, vowelCounter=0, recordingCounter=0;  // looping variables
	float cValuesForVowelA[6*f][p+1] = {{0}}, cValuesForVowelE[6*f][p+1] = {{0}}, cValuesForVowelI[6*f][p+1] = {{0}};
	float cValuesForVowelO[6*f][p+1] = {{0}}, cValuesForVowelU[6*f][p+1] = {{0}};

	int rowCountForA = 0, rowCountForE = 0, rowCountForI = 0,rowCountForO = 0, rowCountForU = 0;

	for (vowelCounter = 0; vowelCounter < 5; vowelCounter++) {
		for (recordingCounter = 1; recordingCounter <= 5; recordingCounter++) {

			int totalSamplesCount = 1, maxSampleValue = 0;
			float sampleSum = 0.0, dcShift = 0.0,normalisedSampleArr[28000], energy[100];
			float stableFrameDataArray[f][N] = {{0}};
			char vowelFilepath[1000], vowelFilename[30];
			
			FILE *vowelSampleFilePtr;
			FILE *vowelSampleFilePtr1;

			strcpy(vowelFilepath, vowelRecordingpath);
			sprintf(vowelFilename, "%c/%s_%c_%d.txt", vowel[vowelCounter], recordingRollNumber, vowel[vowelCounter], recordingCounter);
			strcat(vowelFilepath, vowelFilename);

			fopen_s(&vowelSampleFilePtr, vowelFilepath , "r");
			fopen_s(&vowelSampleFilePtr1, vowelFilepath , "r");
			printf("\n The filePath taken as reference is  %s \n", vowelFilepath);

			for (headerCounter = 0; headerCounter <= 5; headerCounter ++) {
					fgets(headerLine, headerLineSize, vowelSampleFilePtr);
					fgets(headerLine, headerLineSize, vowelSampleFilePtr1);
			}

			// calculating DC shift
			while (fgets(sampleLine, sampleLineSize, vowelSampleFilePtr)) 
			{
				int sampleValue = atoi(sampleLine);
				sampleSum += sampleValue;
				if (maxSampleValue < abs(sampleValue)) 
				{
					maxSampleValue = abs(sampleValue);
				}
				totalSamplesCount++;
			}

			dcShift = sampleSum/totalSamplesCount;

			printf("\n The dc shift  of given sample file %s is %f \n",vowelFilepath, dcShift);
			printf("\n The number of sample count in  file %s is %d \n",vowelFilepath, totalSamplesCount);
	
			// Normalisation of values to +5000 to -5000
			int normalisedValueCounter = 0, maxEnergyIndex = 0;
			float energySum = 0;
			while (fgets(sampleLine, sampleLineSize, vowelSampleFilePtr1)) 
			{
				int sampleValue = atoi(sampleLine);
				float normalisedVAlue = ((float)(sampleValue-dcShift)*5000) / maxSampleValue ;
				normalisedSampleArr[normalisedValueCounter++] = normalisedVAlue;
				if (normalisedValueCounter >= 28000) {
					break;
				}
			}
				
			// converting normalised value 1D Array to 2D Array
			float vowelSampleDataArray[100][N+1] = {{0}};
			int frameNumber = 0;
			int sampleNumber = 0;
			for (i = 0 ; i < normalisedValueCounter; i++) {
				
				if (sampleNumber >= N) {
					frameNumber++;
					sampleNumber=0;
				}
				
				vowelSampleDataArray[frameNumber][sampleNumber++] = normalisedSampleArr[i++];
			}
		
			// Finding energy of each frame
			for(int i=0; i<frameNumber; i++){
				double sum = 0;
				for(int j=0; j<N; ++j){
					sum += (vowelSampleDataArray[i][j] * vowelSampleDataArray[i][j]);
				}
				energy[i] = sum/(N * 1.0);
			}

			// finding the maximum energy frame
			double maxEnergy = 0.0;
			for(int i=0; i<frameNumber; i++){
				if(energy[i] > maxEnergy) {
					maxEnergy = energy[i];
					maxEnergyIndex = i;
				}
			}

			// finding the stable 5 frames
			frameNumber = 0;
			for(int i=maxEnergyIndex-2; i<=maxEnergyIndex+2; i++){
				for(int j=0; j<N; ++j){
					stableFrameDataArray[frameNumber][j] = vowelSampleDataArray[i][j];
				}
				frameNumber++;
			}

			applyHamming(stableFrameDataArray);
	
			// for each frame calculate Ri's, Ai's, Ci's
			for (i = 0; i < f; i++) {

				float refDataArr[N]= {{0}};
				float aArr[p+1] = {{0}};
				float cArr[p+1]= {{0}};
				float rArr[p+1]= {{0}};

				for (j= 0; j< N; j++) {
					refDataArr[j] = stableFrameDataArray[i][j];
				}
				calculateR(refDataArr, rArr,false);
				calculateA(rArr, aArr,false);
				calculateC(cArr, aArr, rArr[0],false);

				if (vowelCounter == 0) {
					for(i=1; i<=p; i++) {
						cValuesForVowelA[rowCountForA][i] = cArr[i];
					}
					rowCountForA++;
				} 
				else if (vowelCounter == 1) {
					for(i=1; i<=p; i++) {
						cValuesForVowelE[rowCountForE][i] = cArr[i];
					}
					rowCountForE++;
				}
				else if (vowelCounter == 2) {
					for(i=1; i<=p; i++) {
						cValuesForVowelI[rowCountForI][i] = cArr[i];
					}
					rowCountForI++;
				}
				else if (vowelCounter == 3) {
					for(i=1; i<=p; i++) {
						cValuesForVowelO[rowCountForO][i] = cArr[i];
					}
					rowCountForO++;
				}
				else if (vowelCounter == 4) {
					for(i=1; i<=p; i++) {
						cValuesForVowelU[rowCountForU][i] = cArr[i];
					}
					rowCountForU++;
				}
			}

		}
	}

	// take the average value of Ci's for the 5 recording of speech
	float cAvgValuesForVowelA[f+1][p+1] = {{0}}, cAvgValuesForVowelE[f+1][p+1] = {{0}}, cAvgValuesForVowelI[f+1][p+1] = {{0}};
	float cAvgValuesForVowelO[f+1][p+1] = {{0}}, cAvgValuesForVowelU[f+1][p+1] = {{0}};

	findAverage(cAvgValuesForVowelA, cValuesForVowelA);
	findAverage(cAvgValuesForVowelE, cValuesForVowelE);
	findAverage(cAvgValuesForVowelI, cValuesForVowelI);
	findAverage(cAvgValuesForVowelO, cValuesForVowelO);
	findAverage(cAvgValuesForVowelU, cValuesForVowelU);

	//apply raise sin Window
	for(int i=0; i<f; ++i){
		for(int j=1; j<=p; ++j){
			cAvgValuesForVowelA[i][j] *= raisedSin(j);
			cAvgValuesForVowelE[i][j] *= raisedSin(j);
			cAvgValuesForVowelI[i][j] *= raisedSin(j);
			cAvgValuesForVowelO[i][j] *= raisedSin(j);
			cAvgValuesForVowelU[i][j] *= raisedSin(j);
		}
	}

	 /*****************************************************************
						 part 3 - Vowel Recognition
	*******************************************************************/


	for (vowelCounter = 0; vowelCounter < 5; vowelCounter++) {
		for (recordingCounter = 6; recordingCounter <= 20 ; recordingCounter++) {

			float cAvgValuesForVowelRecognition[f+1][p+1] = {{0}};
			float cValuesForVowelRecognition[6*f][p+1] = {{0}};
			int rowCountForVowel = 0;
			int totalSamplesCount = 1, maxSampleValue = 0;
			float sampleSum = 0.0, dcShift = 0.0,normalisedSampleArr[28000], energy[100];
			float stableFrameDataArray[f+1][N];
			char vowelFilepath[1000], vowelFilename[30];
			
			FILE *vowelSampleFilePtr;
			FILE *vowelSampleFilePtr1;

			strcpy(vowelFilepath, vowelRecordingpath);
			sprintf(vowelFilename, "%c/%s_%c_%d.txt", vowel[vowelCounter], recordingRollNumber, vowel[vowelCounter], recordingCounter);
			strcat(vowelFilepath, vowelFilename);

			fopen_s(&vowelSampleFilePtr, vowelFilepath , "r");
			fopen_s(&vowelSampleFilePtr1, vowelFilepath , "r");
			
			for (headerCounter = 0; headerCounter <= 5; headerCounter ++) {
					fgets(headerLine, headerLineSize, vowelSampleFilePtr);
					fgets(headerLine, headerLineSize, vowelSampleFilePtr1);
			}


			// calculating DC shift
			while (fgets(sampleLine, sampleLineSize, vowelSampleFilePtr)) 
			{
				int sampleValue = atoi(sampleLine);
				sampleSum += sampleValue;
				if (maxSampleValue < abs(sampleValue)) 
				{
					maxSampleValue = abs(sampleValue);
				}
				totalSamplesCount++;
			}

			dcShift = sampleSum/totalSamplesCount;
			
			printf("\n The dc shift  of given sample file %s is %f \n",vowelFilepath, dcShift);
			printf("\n The number of sample count in  file %s is %d \n",vowelFilepath, totalSamplesCount);
	
			// Normalisation of values to +5000 to -5000
			int normalisedValueCounter = 0, maxEnergyIndex = 0;
			float energySum = 0;
			while (fgets(sampleLine, sampleLineSize, vowelSampleFilePtr1)) 
			{
				int sampleValue = atoi(sampleLine);
				float normalisedVAlue = ((float)(sampleValue-dcShift)*5000) / maxSampleValue ;
				normalisedSampleArr[normalisedValueCounter++] = normalisedVAlue;
				if (normalisedValueCounter >= 28000) {
					break;
				}
			}
				
			
			float vowelSampleDataArray[100][N+1] = {{0}};
			int frameNumber = 0;
			int sampleNumber = 0;
			for (i = 0 ; i < normalisedValueCounter; i++) {
				
				if (sampleNumber >= N) {
					frameNumber++;
					sampleNumber=0;
				}
				
				vowelSampleDataArray[frameNumber][sampleNumber++] = normalisedSampleArr[i++];
			}
		 
			for(int i=0; i<frameNumber; i++){
				double sum = 0;
				for(int j=0; j<N; ++j){
					sum += (vowelSampleDataArray[i][j] * vowelSampleDataArray[i][j]);
				}
				energy[i] = sum/(N * 1.0);
			}

			
			double maxEnergy = 0.0;
			for(int i=0; i<frameNumber; i++){
				if(energy[i] > maxEnergy) {
					maxEnergy = energy[i];
					maxEnergyIndex = i;
				}
			}


			frameNumber = 0;
			for(int i=maxEnergyIndex-2; i<=maxEnergyIndex+2; i++){
				for(int j=0; j<N; ++j){
					stableFrameDataArray[frameNumber][j] = vowelSampleDataArray[i][j];
				}
				frameNumber++;
			}

			applyHamming(stableFrameDataArray);
			for (i = 0; i < f; i++) {

				float refDataArr[N]= {{0}};
				float aArr[p+1] = {{0}};
				float cArr[p+1]= {{0}};
				float rArr[p+1]= {{0}};

				for (j= 0; j< N; j++) {
					refDataArr[j] = vowelSampleDataArray[i][j];
				}
				calculateR(refDataArr, rArr,false);
				calculateA(rArr, aArr,false);
				calculateC(cArr, aArr, rArr[0],false);

				for(i=1; i<=p; i++) {
					cValuesForVowelRecognition[rowCountForVowel][i] = cArr[i];
				}
				rowCountForVowel++;
			}
			
			for(int i=0; i<f; ++i){
				for(int j=1; j<=p; ++j){
					cValuesForVowelRecognition[i][j] *= raisedSin(j);
				}
			}


			double smallestDistance = 9999999.0;
			char vowelRecognised;

			double distanceWithA = computeTokhuraDistance(cValuesForVowelRecognition, cAvgValuesForVowelA, weight);
			if (smallestDistance > distanceWithA) {
				smallestDistance = distanceWithA;
				vowelRecognised = 'A';
			}
			double distanceWithE = computeTokhuraDistance(cValuesForVowelRecognition, cAvgValuesForVowelE, weight);
			if (smallestDistance > distanceWithE) {
				smallestDistance = distanceWithE;
				vowelRecognised = 'E';
			}
			double distanceWithI = computeTokhuraDistance(cValuesForVowelRecognition, cAvgValuesForVowelI, weight);
			if (smallestDistance > distanceWithI) {
				smallestDistance = distanceWithI;
				vowelRecognised = 'I';
			}
			double distanceWithO = computeTokhuraDistance(cValuesForVowelRecognition, cAvgValuesForVowelO, weight);
			if (smallestDistance > distanceWithO) {
				smallestDistance = distanceWithO;
				vowelRecognised = 'O';
			}
			double distanceWithU = computeTokhuraDistance(cValuesForVowelRecognition, cAvgValuesForVowelU, weight);
			if (smallestDistance > distanceWithU) {
				smallestDistance = distanceWithU;
				vowelRecognised = 'U';
			}
			printf("Tokhura distanceWithA = %f, Euclidean distanceWithA = %f \n",distanceWithA, computeEuclideanDistance(cValuesForVowelRecognition, cValuesForVowelA));
			printf("Tokhura distanceWithE = %f, Euclidean distanceWithE = %f  \n",distanceWithE, computeEuclideanDistance(cValuesForVowelRecognition, cValuesForVowelE));
			printf("Tokhura distanceWithI = %f, Euclidean distanceWithI = %f  \n",distanceWithI, computeEuclideanDistance(cValuesForVowelRecognition, cValuesForVowelI));
			printf("Tokhura distanceWithO = %f, Euclidean distanceWithO = %f  \n",distanceWithO, computeEuclideanDistance(cValuesForVowelRecognition, cValuesForVowelO));
			printf("Tokhura distanceWithU = %f, Euclidean distanceWithU = %f  \n",distanceWithU,computeEuclideanDistance(cValuesForVowelRecognition, cValuesForVowelU));
			
			printf("Smallest distance achieved is = %f \n",smallestDistance);
			printf("The vowel recognised is  %c \n", vowelRecognised);
			vowelFilename[0] = '\0';
		}
	}
	return 0;
}



	


