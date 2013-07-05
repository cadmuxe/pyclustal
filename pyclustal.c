#include <Python.h>
#include <stdio.h>
#include "structmember.h"
#include "clustal-omega.h"

static PyObject *AlignError;
static PyObject * pyclustal_getAlign(PyObject *self, PyObject *args){
    const char *filepath;/*{{{*/
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
    if(ReadSequences(prMSeq, filepath, SEQTYPE_UNKNOWN, SQFILE_FASTA, 0, INT_MAX, INT_MAX) == -1){
        PyErr_SetString(PyExc_IOError, "Failed to open the File");
        return NULL;
    }
    if(Align(prMSeq, NULL, &rAlnOpts) == -1){
        PyErr_SetString(AlignError, "Failed to align the sequences");
        return NULL;
    }
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
}/*}}}*/

typedef struct {
    PyObject_HEAD
    opts_t alnOpt; 
} AlignOptions;

static void AlignOptions_dealloc(AlignOptions *self){
    int i = 0;/*{{{*/
    PyMem_Free(self->alnOpt.pcDistmatInfile);
    PyMem_Free(self->alnOpt.pcDistmatOutfile);
    PyMem_Free(self->alnOpt.pcClustfile);
    PyMem_Free(self->alnOpt.pcGuidetreeOutfile);
    PyMem_Free(self->alnOpt.pcGuidetreeInfile);
    for(i = 0; i < self->alnOpt.iHMMInputFiles; i++){
        PyMem_Free(self->alnOpt.ppcHMMInput[i]);
    }
    self->ob_type->tp_free((PyObject *)self);
}/*}}}*/

static PyObject * AlignOptions_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
    AlignOptions *self;/*{{{*/
    self = (AlignOptions *)type->tp_alloc(type, 0);
    if(self != NULL){
        self->alnOpt.bAutoOptions = FALSE;
        self->alnOpt.pcDistmatInfile = NULL;
        self->alnOpt.pcDistmatOutfile = NULL;
        self->alnOpt.iClustersizes = 100;
        self->alnOpt.pcClustfile = NULL;
        self->alnOpt.iClusteringType = 1;
        self->alnOpt.iPairDistType = 1;
        self->alnOpt.bUseMbed = TRUE;
        self->alnOpt.bUseMbedForIteration = TRUE;
        self->alnOpt.pcGuidetreeOutfile = NULL;
        self->alnOpt.pcGuidetreeInfile = NULL;
        self->alnOpt.bPercID = FALSE;
        self->alnOpt.ppcHMMInput = NULL;
        self->alnOpt.iHMMInputFiles = 0;
        self->alnOpt.iNumIterations = 0;
        self->alnOpt.bIterationsAuto = FALSE;
        self->alnOpt.iMaxGuidetreeIterations = INT_MAX;
        self->alnOpt.iMaxHMMIterations = INT_MAX;
        /* 2048 default|give 2GB to MAC algorithm */
        self->alnOpt.rHhalignPara.iMacRamMB = 2048;
        /* protein mode unless we say otherwise */
        self->alnOpt.rHhalignPara.bIsDna = false;  
        self->alnOpt.rHhalignPara.bIsRna = false;    
    }
    return (PyObject *)self;
}/*}}}*/

static PyMemberDef AlignOptions_members[]={
    {"autoOptions", T_INT, /*{{{*/
        offsetof(AlignOptions,alnOpt) + offsetof(opts_t, bAutoOptions),
        0,
        "value:(0,1) Clustal (know what) is good for you"},
    {"clusteringType", T_INT,
        offsetof(AlignOptions,alnOpt) + offsetof(opts_t, iClusteringType),
        0,
        "clustering type"},
    {"clusterSizes", T_INT,
        offsetof(AlignOptions,alnOpt) + offsetof(opts_t, iClustersizes),
        0,
        "number of sequences in cluster"},
    {"pairDistType", T_INT,
        offsetof(AlignOptions,alnOpt) + offsetof(opts_t, iPairDistType),
        0,
        "pairwise distance method"},
    {"useMbed", T_INT,
        offsetof(AlignOptions,alnOpt) + offsetof(opts_t, bUseMbed),
        0,
        "value:(0,1) use mbed-link clustering"},
    {"useMbedForIteration", T_INT,
        offsetof(AlignOptions,alnOpt) + offsetof(opts_t, bUseMbedForIteration),
        0,
        "value:(0,1) use mbed-like clustering also during iteration"},
    {"useKimra", T_INT,
        offsetof(AlignOptions,alnOpt) + offsetof(opts_t, bUseKimura),
        0,
        "value:(0,1) use Kimura corrected distance"},
    {"percID", T_INT,
        offsetof(AlignOptions,alnOpt) + offsetof(opts_t, bPercID),
        0,
        "value:(0,1) print percentage identity"},
    {"numIterations", T_INT,
        offsetof(AlignOptions,alnOpt) + offsetof(opts_t, iNumIterations),
        0,
        "number of iterations"},
    {"iterationsAuto", T_INT,
        offsetof(AlignOptions,alnOpt) + offsetof(opts_t, bIterationsAuto),
        0,
        "value:(0,1) determine number of iterations automatically"},
    {"maxHMMIterations", T_INT,
        offsetof(AlignOptions,alnOpt) + offsetof(opts_t, iMaxHMMIterations),
        0,
        "maximum number of hmm iterations"},
    {"maxGuidetreeIterations", T_INT,
        offsetof(AlignOptions,alnOpt) + offsetof(opts_t, iMaxGuidetreeIterations),
        0,
        "max number of guidetree iterations"},
    {"MacRamMB", T_INT,
        offsetof(AlignOptions,alnOpt) + offsetof(opts_t, rHhalignPara) + 
            offsetof(hhalign_para, iMacRamMB),
        0,
        "dedicated amount of RAM for Maximum Accuracy (in MB)"},
    {"isDna", T_INT,
        offsetof(AlignOptions,alnOpt) + offsetof(opts_t, rHhalignPara) +
            offsetof(hhalign_para, bIsDna),
        0,
        "value:(0,1) indicates we're in nucleotide mode"},
    {"isRna", T_INT,
        offsetof(AlignOptions,alnOpt) + offsetof(opts_t, rHhalignPara) +
            offsetof(hhalign_para, bIsRna),
        0,
        "value:(0,1) indicates we're in nucleotide mode"},
    {NULL}
};/*}}}*/

static PyMethodDef AlignOptions_methods[]={
    {NULL}
};
static PyTypeObject AlignOptionsType={
    PyObject_HEAD_INIT(NULL)
    0,                          /*ob_size*/ /*{{{*/
    "pyclustal.AlignOptions",              /*tp_name*/
    sizeof(AlignOptions),              /*tp_basicsize*/ 
    0,                          /*tp_itemsize*/ 
    (destructor)AlignOptions_dealloc,  /*tp_dealloc*/ 
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    0,                          /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash*/
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,         /*tp_flags*/
    "Alignment options",            /*tp_doc*/
    0,                          /*tp_traverse*/
    0,                          /*tp_clear*/
    0,                          /*tp_richcompare*/
    0,                          /*tp_weaklistoffset*/
    0,                          /*tp_iter*/
    0,                          /*tp_iternext*/
    AlignOptions_methods,              /*tp_methods*/
    AlignOptions_members,              /*tp_members*/
    0,                          /*tp_getset*/
    0,                          /*tp_base*/
    0,                          /*tp_dict*/
    0,                          /*tp_descr_get*/
    0,                          /*tp_descr_set*/
    0,                          /*tp_dictoffset*/
    0,       /*tp_init*/
    0,                          /*tp_alloc*/
    AlignOptions_new,                  /*tp_new*/
};/*}}}*/

static PyMethodDef pyclustal_methods[] = {
    {"getAlign", (PyCFunction) pyclustal_getAlign, METH_VARARGS,
    "getAlign(filename) -> result \nGet Multiple sequence alignment. "},
    {NULL, NULL, 0, NULL}
};

static int uti_convert_string(char **ppchar, PyObject *value){
    if (!PyString_Check(value)){
        PyErr_SetString(PyExc_TypeError, "The attribute must be string");
    }
    if (value == NULL){
        PyErr_SetString(PyExc_TypeError, "The attribute can not be None");
        return -1;
    }
    if( PyString_Size(value) == 0){
        PyMem_Free(*ppchar);
        (*ppchar) = NULL;
    }
    else{
        ;
    }
}

static int AlignOptions_setDistmatInfile(AlignOptions *self, PyObject *value, void *closure){
    
}

PyMODINIT_FUNC initpyclustal(void){
    PyObject *m;
    if(PyType_Ready(&AlignOptionsType) < 0)
        return;
    m = Py_InitModule3("pyclustal", pyclustal_methods,
            "A Clustal-Omega wrap.");
    if (m == NULL)
        return;
    AlignError = PyErr_NewException("pyclustal.AlignError", NULL, NULL);
    Py_INCREF(AlignError);
    PyModule_AddObject(m, "AlignError", AlignError);
    PyModule_AddObject(m, "AlignOptions", (PyObject *)&AlignOptionsType);
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
