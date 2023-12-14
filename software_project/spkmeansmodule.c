#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <math.h>

#include "spkmeans.h"




/* --------------------API FROM HERE----------------------------------*/

static PyObject* new_fit_Kmeans(PyObject *self, PyObject *args)
{
    PyObject *list, *item, *listToPython;
    Py_ssize_t i, n, listPythonLength;
    
    long numberOfPoints, pointSize, k,tmp_long;
    int check;
    int goal;
    double *listFromKmean;
    
    if(!PyArg_ParseTuple(args, "Ollli", &list, &numberOfPoints, &pointSize, &k, &goal)) {
        return NULL;
    }
    /* Is it a list? */
    if (!PyList_Check(list))
        Py_RETURN_NONE;
    
    /*  size of list from Python */
    n = PyList_Size(list);
    
    double *listToKmean = malloc(sizeof listToKmean * n);
    if(listToKmean==NULL){
        Py_RETURN_NONE;
    }
    
    for (i = 0; i < n; i++) {
        item = PyList_GetItem(list, i);
        listToKmean[i] =  PyFloat_AsDouble(item);

    }
    check =  kmean_new_fit(&listFromKmean ,listToKmean, &numberOfPoints, &pointSize, &k,  goal);
    
    free(listToKmean);
    
    if(check==1) {
        Py_RETURN_NONE;
    }
    
    else{
        tmp_long = (pointSize * numberOfPoints) + 2;
        
        
        listPythonLength = PyLong_AsSsize_t(PyLong_FromLong(tmp_long));
        
        listToPython = PyList_New(listPythonLength);
        
        for(i=0; i < listPythonLength-2; i++){
            PyList_SetItem(listToPython, i ,PyFloat_FromDouble(listFromKmean[i]));
        }
        
        PyList_SetItem(listToPython, listPythonLength - 2 ,PyFloat_FromDouble( (double)pointSize ) );
        
        listPythonLength = listPythonLength -1;
        PyList_SetItem(listToPython, (listPythonLength) ,PyFloat_FromDouble( (double)k ) );
        
        free(listFromKmean);
    }
    return  listToPython;
}





/*-------------------*****  2nd stage function, no need to change ****-----------------------*/


static PyObject* fit_Kmeans(PyObject *self, PyObject *args)
{
    PyObject *list, *item, *listToPython;
    Py_ssize_t i, n, listPythonLength;
    
    double *plist;
    long numberOfPoints, pointSize, k;
    int maxIter;
    double epsilon;
    int check;
    if(!PyArg_ParseTuple(args, "Olllidn", &list, &numberOfPoints, &pointSize, &k, &maxIter, &epsilon, &listPythonLength)) {
        return NULL;
    }
    /* Is it a list? */
    if (!PyList_Check(list))
        Py_RETURN_NONE;;
    /* Get the size of it and build the output list */
    n = PyList_Size(list);
    
    double *Clist = malloc(sizeof Clist * n);
    if(Clist==NULL){
        Py_RETURN_NONE;
    }
    
    for (i = 0; i < n; i++) {
        item = PyList_GetItem(list, i);
        Clist[i] =  PyFloat_AsDouble(item);

    }
    check =  kmean_fit(Clist, numberOfPoints, pointSize, k, maxIter, epsilon, &plist);
    
    free(Clist);
 
    if(check==1){
        Py_RETURN_NONE;
    }
    
    else{
        listToPython = PyList_New(listPythonLength);
        for(i=0; i<(listPythonLength); i++){
            PyList_SetItem(listToPython, i ,PyFloat_FromDouble(plist[i]));
            
        }
        free(plist);
    }
    return  listToPython;
}



/*--------------*****  end of API functions   ****------------------*/


static PyMethodDef capiMethods[] = {
    {"fit",                   /* the Python method name that will be used */
      (PyCFunction) fit_Kmeans, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parametersaccepted for this function */
      PyDoc_STR("the kmeans alg")}, /*  The docstring for the function */
    
    {"new_fit", (PyCFunction) new_fit_Kmeans, METH_VARARGS, PyDoc_STR("the kmeans alg")},
    {NULL, NULL, 0, NULL}     /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};


/* This initiates the module using the above definitions. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,"mykmeanssp", NULL, -1, capiMethods
};


/*
 * The PyModuleDef structure, in turn, must be passed to the interpreter in the module’s initialization function.
 * The initialization function must be named PyInit_name(), where name is the name of the module and should match
 * what we wrote in struct PyModuleDef.
 * This should be the only non-static item defined in the module file
 */
PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}




