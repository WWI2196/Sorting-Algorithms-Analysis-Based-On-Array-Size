#include <iostream>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <cstring>
#include <string>

using namespace std;
using namespace std::chrono;
// added a comment

const int NAMECOLUMNWIDTH = 25,NCOLUMNWIDTH = 10;// Column widths for the table

const int NVALUETYPES = 6, NUMBEROFRUNS = 10; // Number of values for n and number of runs
const int NVALUES[NVALUETYPES] = { 100, 200, 400, 800, 1000, 2000 }; // Array of values for n
const string SORTINGNAMES[] = { "Selection Sort","Insertion Sort","Merge Sort","Quick Sort","Counting Sort2" }; // Array of sorting algorithm names

void printTableHeader();
void swap(int* a, int* b);
void generateRandomArray(int* numberArray, int size);
void printAverageComparisonsTable(const string sortTypeName,int* NVALUES, long long comparisons[], int NVALUETYPES);
void printAverageRunTimeTable(const string sortTypeName, int* NVALUES, double runTimes[], int NVALUETYPES);
void selectionSort(int numberArray[], int size, double& runTime, long long& comparisons);
void insertionSort(int numberArray[], int size, double& runTime, long long& comparisons);
void merge(int arr[], int l, int m, int r, long long& comparisons);
void mergeSort(int numberArray[], int l, int r, long long& comparisons);
int partition(int numberArray[], int low, int high, long long& comparisons);
void quickSort(int numberArray[], int low, int high, long long& comparisons);
void countingSort_R(int numberArray[], int size, double& runTime, long long& comparisons);
void sortSelection(int* numberArray, int arraySize, int sortType, long long& comparisons, double& runTime);
void sortingAverageValues(const int* NVALUES, const int NVALUETYPES, int NUMBEROFRUNS);

void printTableHeader() {// Print the table header

    cout << '|' << string(91, '-') << '|' << endl;
    
    cout << '|' << left << setw(NAMECOLUMNWIDTH) << "Sorting Algorithm" << '|' << setw(NCOLUMNWIDTH) << "n=100" << '|' << setw(NCOLUMNWIDTH) << "n=200" 
        << '|' << setw(NCOLUMNWIDTH) << "n=400"<< '|' << setw(NCOLUMNWIDTH) << "n=800" << '|' << setw(NCOLUMNWIDTH) << "n=1000" << '|' << setw(NCOLUMNWIDTH) << "n=2000" << '|' << endl;

    cout << '|' << string(NAMECOLUMNWIDTH, '-') << '|' << string(NCOLUMNWIDTH, '-') << '|' << string(NCOLUMNWIDTH, '-') << '|' << string(NCOLUMNWIDTH, '-') 
        << '|' << string(NCOLUMNWIDTH, '-') << '|' << string(NCOLUMNWIDTH, '-') << '|' << string(NCOLUMNWIDTH, '-') << '|' << endl;
}

void printAverageComparisonsTable(string sortTypeName,const int* NVALUES, long long comparisons[],const int NVALUETYPES) { // print the average number of comparisons table
    cout << '|' << left << setw(NAMECOLUMNWIDTH) << sortTypeName;

    for (int i = 0; i < NVALUETYPES; i++) {
        cout << '|' << setw(NCOLUMNWIDTH) << comparisons[i]; // for each value of n print the average number of comparisons
    }
    cout << '|' << endl;

}

void printAverageRunTimeTable(string sortTypeName, const int* NVALUES, double runTimes[],const int NVALUETYPES) { // print the average running time table
    cout << '|' << left << setw(NAMECOLUMNWIDTH) << sortTypeName;

    for (int i = 0; i < NVALUETYPES; i++) {
        cout << '|' << setw(NCOLUMNWIDTH) << runTimes[i]; // for each value of n print the average running time
    }
    cout << '|' << endl;
}

void swap(int* a, int* b) { // swap two number positions
    int temp = *a;
    *a = *b;
    *b = temp;
}

void generateRandomArray(int* numberArray, int size) { // generate random numbers for the array
    srand(time(0)); // seed the random number generator to get different numbers each time

    for (int i = 0; i < size; i++) {
        numberArray[i] = (rand() % (size * 10)) + 1; // generate random numbers between 1 and 10*size and add to the number array
    }
}

void selectionSort(int numberArray[], int size, double& runTime, long long& comparisons) { //function to sort the array using selection sort
    auto start = high_resolution_clock::now();

    for (int i = 0; i < size - 1; i++) {
        for (int y = i + 1; y < size; y++) {
            comparisons++;
            if (numberArray[i] > numberArray[y]) swap(numberArray[i], numberArray[y]);
        }
    }

    auto end = high_resolution_clock::now();
    runTime = duration_cast<chrono::duration<double, milli>>(end - start).count();
}

void insertionSort(int numberArray[], int size, double& runTime, long long& comparisons) { //function to sort the array using insertion sort
    auto start = high_resolution_clock::now();

    for (int i = 1; i < size; i++) {
        int j = i;
        comparisons++;

        while (j > 0 && numberArray[j] < numberArray[j - 1]) {
            swap(&numberArray[j], &numberArray[j - 1]);
            j--;
            comparisons++;
        }
    }

    auto end = high_resolution_clock::now();
    runTime = duration_cast<chrono::duration<double, milli>>(end - start).count();
}

void merge(int numberArray[], int l, int m, int r, long long& comparisons) { // Merge function for merge sort
    int i, j;
    int n1 = m - l + 1;
    int n2 = r - m;

    int* L = new int[n1];  // Dynamically allocate memory for L
    int* R = new int[n2];  // Dynamically allocate memory for R

    for (i = 0; i < n1; i++) {
        L[i] = numberArray[l + i];
        comparisons++;
    }

    for (j = 0; j < n2; j++) {
        R[j] = numberArray[m + 1 + j];
        comparisons++;
    }

    i = 0;
    j = 0;

    while (i < n1 && j < n2) {
        comparisons++;
        if (L[i] <= R[j]) {
            numberArray[l] = L[i];
            i++;
        }
        else {
            numberArray[l] = R[j];
            j++;
        }
        l++;
    }

    while (i < n1) {
        numberArray[l] = L[i];
        i++;
        l++;
        comparisons++;
    }

    while (j < n2) {
        numberArray[l] = R[j];
        j++;
        l++;
        comparisons++;
    }

    delete[] L;  // Deallocate memory for L
    delete[] R;  // Deallocate memory for R
}

void mergeSort(int numberArray[], int l, int r, long long& comparisons) { // function for merge sort
    if (l < r) {
        comparisons++;
        int m = l + (r - l) / 2;

        mergeSort(numberArray, l, m, comparisons);
        mergeSort(numberArray, m + 1, r, comparisons);

        merge(numberArray, l, m, r, comparisons);
    }
}

int partition(int numberArray[], int low, int high, long long& comparisons) { // function for partitioning the array for quick sort
    int pivot = numberArray[high];
    int i = (low - 1);

    for (int j = low; j < high; j++) {
        comparisons++;
        if (numberArray[j] <= pivot) {
            i++;
            swap(&numberArray[i], &numberArray[j]);
        }
    }

    comparisons++;
    swap(&numberArray[i + 1], &numberArray[high]);

    return (i + 1);
}

void quickSort(int numberArray[], int low, int high, long long& comparisons) { // function for quick sort
    if (low < high) {
        comparisons++;
        int pivot = partition(numberArray, low, high, comparisons);
        quickSort(numberArray, low, pivot - 1, comparisons);
        quickSort(numberArray, pivot + 1, high, comparisons);
    }
}

void countingSort_R(int numberArray[], int size, double& runTime, long long& comparisons) { // function for counting sort
    auto start = high_resolution_clock::now();

    int* index_array = new int[size]();  // Initialize with zero
    int* sorted_array = new int[size];

    for (int i = 0; i < size; i++) {
        comparisons++;
        int count = 0;
        for (int y = 0; y < size; y++) {
            comparisons++;
            if (numberArray[y] < numberArray[i]) {
                count++;
            }
        }
        index_array[i] = count;
    }

    for (int i = 0; i < size; i++) {
        comparisons++;
        sorted_array[index_array[i]] = numberArray[i];
    }

    for (int i = 0; i < size; i++) {
        comparisons++;
        numberArray[i] = sorted_array[i];
    }

    auto end = high_resolution_clock::now();
    runTime = duration_cast<chrono::duration<double, milli>>(end - start).count();

    delete[] index_array;
    delete[] sorted_array;

}

void sortSelection(int* numberArray, int arraySize, int sortType, long long& comparisons, double& runTime) { // function to select the sorting algorithm based on the user's choice
    int* copyArray = new int[arraySize]; // array to store the copy of the original array
    copy(numberArray, numberArray + arraySize, copyArray); // copy the original number array to the copy array

    switch (sortType) {
    case 1:
        selectionSort(copyArray, arraySize, runTime, comparisons);
        break;
    case 2:
        insertionSort(copyArray, arraySize, runTime, comparisons);
        break;
    case 3:
    {
        auto start = high_resolution_clock::now();
        mergeSort(copyArray, 0, arraySize - 1, comparisons);
        auto end = high_resolution_clock::now();
        runTime = duration_cast<chrono::duration<double, milli>>(end - start).count();
        break;
    }
    case 4:
    {
        auto start = high_resolution_clock::now();
        quickSort(copyArray, 0, arraySize - 1, comparisons);
        auto end = high_resolution_clock::now();
        runTime = duration_cast<chrono::duration<double, milli>>(end - start).count();
        break;
    }
    case 5:
        countingSort_R(copyArray, arraySize, runTime, comparisons);
        break;
    default:
        break;
    }

    delete[] copyArray;
}

void sortingAverageValues(const int* NVALUES, const int NVALUETYPES, int NUMBEROFRUNS) {
    long long comparisons[5][6] = { 0 }; // 2D array to store the average number of comparisons for each sorting process and n value
    double runTimes[5][6] = { 0 }; // 2D array to store the average running time for each sorting process and n value

    for (int i = 0; i < NVALUETYPES; i++) { // for each type of n value
        int* numberArray = new int[NVALUES[i]]; // craete an array of given n value

        for (int run = 0; run < NUMBEROFRUNS; run++) { // for each run 
            generateRandomArray(numberArray, NVALUES[i]); // generate random numbers for the array 

            for (int sortType = 1; sortType <= 5; sortType++) { // for each sorting type
                long long* comparison_temp = new long long(0); // temporary variable to store the number of comparisons
                double* runTime_temp = new double(0); // temporary variable to store the running time

                sortSelection(numberArray, NVALUES[i], sortType, *comparison_temp, *runTime_temp); // generate the values for the temporary variables
                comparisons[sortType - 1][i] += *comparison_temp; // add the temporary values to the 2D array for each run 
                runTimes[sortType - 1][i] += *runTime_temp; // add the temporary values to the 2D array for each run

                delete comparison_temp;
                delete runTime_temp;
            }
        }
        delete[] numberArray;
    }

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < NVALUETYPES; j++) {
            comparisons[i][j] /= NUMBEROFRUNS; // calculate the average number of comparisons for each sorting process and n value by dinding the cell value by number of runs
            runTimes[i][j] /= NUMBEROFRUNS; // calculate the average running time for each sorting process and n value by dinding the cell value by number of runs
        }
    }

    cout << "\tAverage number of comparisons for sorting arrays of n integers" << endl;
    printTableHeader(); // print the table header
    // print the table for the average number of comparisons
    for (int i = 0; i < 5; i++) { 
        printAverageComparisonsTable(SORTINGNAMES[i], NVALUES, comparisons[i], NVALUETYPES); // print the average comparisons valus for each sorting types and n values
    }
    cout << endl;

    cout << "\tAverage running time(in ms) for sorting arrays of n integers" << endl;
    printTableHeader(); // print the table header
    // print the table for the average running time
    for (int i = 0; i < 5; i++) {
        printAverageRunTimeTable(SORTINGNAMES[i], NVALUES, runTimes[i], NVALUETYPES); // print the average running time for each sorting types and n values
    }
    cout << endl;
}

int main() {
    
    sortingAverageValues(NVALUES, NVALUETYPES, NUMBEROFRUNS); // call the function to calculate the average number of comparisons and running time for each sorting process and n value

    return 0;
}
