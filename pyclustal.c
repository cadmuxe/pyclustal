#include <Python.h>
#include <stdio.h>
#include "clustal-omega.h"

static PyObject * pyclustal_getAlign(PyObject *self, PyObject *args){
    const char *filepath;
    mseq_t *prMSeq = NULL;
    int iThreads = 1;
    opts_t rAlnOpts;
    int iAux;
    PyObject *result = NULL;
    if (!PyArg_ParseTuple(args, "s", &filepath))
        return NULL;
    
    LogDefaultSetup(&rLog);
    SetDefaultAlnOpts(&rAlnOpts);
    InitClustalOmega(iThreads);
    NewMSeq(&prMSeq);
    ReadSequences(prMSeq, filepath, SEQTYPE_UNKNOWN, SQFILE_FASTA, 0, INT_MAX, INT_MAX);
    Align(prMSeq, NULL, &rAlnOpts);
    result = PyList_New(prMSeq->nseqs);
    for(iAux = 0; iAux < prMSeq->nseqs; iAux++){
        PyList_SetItem(result, iAux, Py_BuildValue("{s:s,s:s,s:s,s:s}",
                    "name",     prMSeq->sqinfo[iAux].name,
                    "desc",     prMSeq->sqinfo[iAux].desc,
                    "origin",   prMSeq->orig_seq[iAux],
                    "seq",      prMSeq->seq[iAux]
                    ));
    }
    FreeMSeq(&prMSeq);
    return result;
}
static PyMethodDef pyclustal_methods[] = {
    {"getAlign", (PyCFunction) pyclustal_getAlign, METH_VARARGS,
    "getAlign(filename) -> result \nGet Multiple sequence alignment. "},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initpyclustal(void){
    Py_InitModule("pyclustal", pyclustal_methods);
}


/*
int main(int argc, char **argv){
    mseq_t *prMSeq = NULL;
    int iThreads = 1;
    opts_t rAlnOpts;
    int iAux;
    char *pfile = "test.fa";
    LogDefaultSetup(&rLog);
    SetDefaultAlnOpts(&rAlnOpts);
    InitClustalOmega(iThreads);
    printf("GO\n");
    NewMSeq(&prMSeq);
    ReadSequences(prMSeq, pfile, SEQTYPE_UNKNOWN, SQFILE_FASTA, 0, INT_MAX, INT_MAX);
    Align(prMSeq, NULL, &rAlnOpts);
    printf("%d\n", prMSeq->nseqs);
    for(iAux = 0; iAux < prMSeq->nseqs; iAux++){
        printf("%d:\t %s\n",iAux, prMSeq->seq[iAux]);
    }
    printf("Begin!\n");
    FreeMSeq(&prMSeq);
    return 0;
}
*/