#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "spkmeans.h"



static PyObject* getDataPoints_capi(PyObject *self, PyObject *args);
static void convert_cMatrix_to_pyMAtrix(PyObject* pyMatrix, double** cMatrix, int n,int d);
static void convert_pyMatrix_to_cMatrix(PyObject* pyMatrix, double** cMatrix, int n, int d);
static PyObject* wam(PyObject *self, PyObject *args);

static PyObject* getDataPoints_capi(PyObject *self, PyObject *args){
    /**
     * @brief for an input of python observations return the T (n X k) matrix as part of the SPH proccess
     * 
     */
   
    PyObject* source_dataPy;
    PyObject* result_dataPy;
    int n,d,k;
    double** sourceDataForC, **newDataPoints;

    if(!PyArg_ParseTuple(args,"iiiO",&n,&d,&k,&source_dataPy)){ /*recieve the python arguments*/
        return NULL;
    }

    result_dataPy = PyList_New(n);
    sourceDataForC = allocationMatrix(n,d);
    convert_pyMatrix_to_cMatrix(source_dataPy,sourceDataForC,n,d); /* convert the matrix to c format*/
    newDataPoints = getDataPoints(sourceDataForC,n,d,&k); /* using the spkmeans process to get tthe new data points*/
    convert_cMatrix_to_pyMAtrix(result_dataPy,newDataPoints,n,k); /*convert the data back to python format*/
    freeMatrix(sourceDataForC,n);
    return result_dataPy;

}

static PyObject* kmeansAlgorithm_capi(PyObject *self, PyObject *args)
{
    /**
     * @brief for an input of observations data (n X k) and K init clusters, perform the kmenas algorithm 
     * 
     */
    PyObject* observations_py;
    PyObject* centroids_py;
    int n,d,k;
    double** observations_c;
    double** centroids;

    if(!PyArg_ParseTuple(args,"iiiOO",&n,&d,&k,&observations_py,&centroids_py)){
        return NULL;
    }

    observations_c = allocationMatrix(n,k);
    centroids = allocationMatrix(k,k);
    convert_pyMatrix_to_cMatrix(observations_py,observations_c,n,d);
    convert_pyMatrix_to_cMatrix(centroids_py,centroids,k,k);
    kmeansAlgorithm(observations_c,centroids,n,d,k);
    freeMatrix(centroids,k);
    freeMatrix(observations_c,n);
    Py_RETURN_NONE;

}

static PyObject* wam(PyObject *self, PyObject *args){
    /**
     * @brief for an input of observations form and print the weighted adjacency matrix
     * 
     */
    PyObject* source_dataPy;
    int n,d;
    double** sourceDataForC;

    if(!PyArg_ParseTuple(args,"iiO",&n,&d,&source_dataPy)){
        return NULL;
    }
    sourceDataForC = allocationMatrix(n,d);
    convert_pyMatrix_to_cMatrix(source_dataPy,sourceDataForC,n,d);
    wamProcess(sourceDataForC, n,d);
    freeMatrix(sourceDataForC,n);
    Py_RETURN_NONE;
}

static PyObject* ddg(PyObject *self, PyObject *args){
    /**
     * @brief for an input of observations form and print the Diagonal Degree matrix
     * 
     */
    PyObject* source_dataPy;
    int n,d;
    double** sourceDataForC;

    if(!PyArg_ParseTuple(args,"iiO",&n,&d,&source_dataPy)){
        return NULL;
    }
    sourceDataForC = allocationMatrix(n,d);
    convert_pyMatrix_to_cMatrix(source_dataPy,sourceDataForC,n,d);
    ddgProcess(sourceDataForC, n,d);
    freeMatrix(sourceDataForC,n);
    Py_RETURN_NONE;
}

static PyObject* jacobi(PyObject *self, PyObject *args){

    /**
     * @brief for an input of observations form and print the the eigenValues and vectors -> the output of jacobi algorithm 
     * 
     */
    PyObject* source_dataPy;
    int n,d;
    double** sourceDataForC;

    if(!PyArg_ParseTuple(args,"iiO",&n,&d,&source_dataPy)){
        return NULL;
    }
    sourceDataForC = allocationMatrix(n,d);
    convert_pyMatrix_to_cMatrix(source_dataPy,sourceDataForC,n,d);
    jacobiProcess(sourceDataForC, n);
    freeMatrix(sourceDataForC,n);
    Py_RETURN_NONE;
}

static PyObject* lnorm(PyObject *self, PyObject *args){

    /**
     * @brief for an input of observations form and print the Normalized Graph Laplacian
     * 
     */
    PyObject* source_dataPy;
    int n,d;
    double** sourceDataForC;

    if(!PyArg_ParseTuple(args,"iiO",&n,&d,&source_dataPy)){
        return NULL;
    }
    sourceDataForC = allocationMatrix(n,d);
    convert_pyMatrix_to_cMatrix(source_dataPy,sourceDataForC,n,d);
    lnormProcess(sourceDataForC, n,d);
    freeMatrix(sourceDataForC,n);
    Py_RETURN_NONE;
}

static void convert_pyMatrix_to_cMatrix(PyObject* pyMatrix, double** cMatrix, int n, int d){
    int i;
    int j;
    double cItem;
    PyObject* py_point;
    double* c_point;

    for (i = 0; i< n; i++){
        py_point = PyList_GET_ITEM(pyMatrix,i);
        c_point = allocationVector(d);
        for (j = 0; j < d; j++){c_point[j] =  PyFloat_AsDouble(PyList_GetItem(py_point,j));}

        for(j=0; j <d ; j++){
            cItem = c_point[j];
            cMatrix[i][j] = cItem;
        }
        free(c_point);  
    }
}


static void convert_cMatrix_to_pyMAtrix(PyObject* pyMatrix, double** cMatrix, int n,int d){

    int i;
    int j;
    PyObject* line;

    for(i = 0; i < n; i++){
        line = PyList_New(d);
        for (j = 0; j < d; j++){
            PyList_SetItem(line,j,PyFloat_FromDouble(cMatrix[i][j]));
        }
        PyList_SetItem(pyMatrix,i,line);
    }

    
}

static PyMethodDef _capiMethods[] = { 
    {"wam", (PyCFunction) wam, METH_VARARGS, PyDoc_STR("form WAM")},
    {"ddg", (PyCFunction) ddg, METH_VARARGS, PyDoc_STR("form DDG")},
    {"lnorm", (PyCFunction) lnorm, METH_VARARGS, PyDoc_STR("form LNORM")},
    {"jacobi", (PyCFunction) jacobi, METH_VARARGS, PyDoc_STR("perform jacobi algorithm - output eignvalues & eignVectors")},
    {"get_new_data", (PyCFunction) getDataPoints_capi, METH_VARARGS, PyDoc_STR("process the input threw the SPK process to get the new datapoionts")},
    {"kmeans", (PyCFunction) kmeansAlgorithm_capi, METH_VARARGS, PyDoc_STR("calculate kmeans from given datapoints")},
    {NULL,NULL,0,NULL}
};

static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeans",
    NULL,
    -1,
    _capiMethods
};

PyMODINIT_FUNC
PyInit_spkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&_moduledef);
    if(!m){
        return NULL;
    }
    return m;
}