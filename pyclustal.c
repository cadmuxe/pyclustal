#include <Python.h>
#include <stdio.h>
#include "structmember.h"
#include "clustal-omega.h"

static PyObject *AlignError;
static PyObject * pyclustal_getAlign(PyObject *self, PyObject *args){
    char *filepath;/*{{{*/
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

static PyObject * pyclustal_Align(PyObject *self, PyObject *args){
    char *filepath;/*{{{*/
    mseq_t *prMSeq = NULL;
    int iThreads = 1;
    int iAux;
    AlignOptions *opt = NULL;
    PyObject *result = NULL;
    if (!PyArg_ParseTuple(args, "sO", &filepath, &opt))
        return NULL;

    printf("OK\n");
    fflush(stdout);

    LogDefaultSetup(&rLog);
    InitClustalOmega(iThreads);
    NewMSeq(&prMSeq);
    if(ReadSequences(prMSeq, filepath, SEQTYPE_UNKNOWN, SQFILE_FASTA, 0, INT_MAX, INT_MAX) == -1){
        PyErr_SetString(PyExc_IOError, "Failed to open the File");
        return NULL;
    }
    if(Align(prMSeq, NULL, &(opt->alnOpt)) == -1){
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

}


static void AlignOptions_dealloc(AlignOptions *self){
    /*{{{*/ 
    int i = 0;
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

static int uti_PyString_to_pchar(char **ppchar, PyObject *value){
    int size = 0;/*{{{*/
    if (!PyString_Check(value)){
        PyErr_SetString(PyExc_TypeError, "The attribute must be string");
        return -1;
    }
    if (value == NULL){
        PyErr_SetString(PyExc_TypeError, "The attribute can not be None");
        return -1;
    }
    size = PyString_Size(value) + 1; // extra one is for sentinel
    if (size == 0){
        PyMem_Free(*ppchar);
        (*ppchar) = NULL;
    }
    else{
        if((*ppchar) == NULL)
            (*ppchar) = PyMem_New(char, size);
        else
            (*ppchar) = PyMem_Resize((*ppchar), char, size);
        memcpy((*ppchar), PyString_AsString(value), size);
    }
    return 0;
}/*}}}*/

static PyObject *uti_pchar_toPyString(char *pchar){
    if(pchar == NULL)
        PyString_FromString("");
    else
        return PyString_FromString(pchar);
}

static int AlignOptions_setDistmatInfile(AlignOptions *self, PyObject *value, void *closure){
    if(uti_PyString_to_pchar(&(self->alnOpt.pcDistmatInfile), value) == -1)
        return -1;
    return 0;
}
static PyObject * AlignOptions_getDistmatInfile(AlignOptions *self, void *closure){
    return uti_pchar_toPyString(self->alnOpt.pcDistmatInfile);
}

static int AlignOptions_setDistmatOutfile(AlignOptions *self, PyObject *value, void *closure){
    if(uti_PyString_to_pchar(&(self->alnOpt.pcDistmatOutfile), value) == -1)
        return -1;
    return 0;
}
static PyObject * AlignOptions_getDistmatOutfile(AlignOptions *self, void *closure){
    return uti_pchar_toPyString(self->alnOpt.pcDistmatOutfile);
}

static int AlignOptions_setClustfile(AlignOptions *self, PyObject *value, void *closure){
    if(uti_PyString_to_pchar(&(self->alnOpt.pcClustfile), value) == -1)
        return -1;
    return 0;
}
static PyObject * AlignOptions_getClustfile(AlignOptions *self, void *closure){
    return uti_pchar_toPyString(self->alnOpt.pcClustfile);
}

static int AlignOptions_setGuidetreeOutfile(AlignOptions *self, PyObject *value, void *closure){
    if(uti_PyString_to_pchar(&(self->alnOpt.pcGuidetreeOutfile), value) == -1)
        return -1;
    return 0;
}
static PyObject * AlignOptions_getGuidetreeOutfile(AlignOptions *self, void *closure){
    return uti_pchar_toPyString(self->alnOpt.pcGuidetreeOutfile);
}

static int AlignOptions_setGuidetreeInfile(AlignOptions *self, PyObject *value, void *closure){
    if(uti_PyString_to_pchar(&(self->alnOpt.pcGuidetreeInfile), value) == -1)
        return -1;
    return 0;
}
static PyObject * AlignOptions_getGuidetreeInfile(AlignOptions *self, void *closure){
    return uti_pchar_toPyString(self->alnOpt.pcGuidetreeInfile);
}

static int AlignOptions_setHMMInput(AlignOptions *self, PyObject *value, void *closure){
    int i,size;
	size = self->alnOpt.iHMMInputFiles;
	// Free all the space
	for(i = 0; i < size; i++){
		PyMem_Free(self->alnOpt.ppcHMMInput[i]);
	}
	PyMem_Free(self->alnOpt.ppcHMMInput);

	size = PyTuple_Size(value);
    printf("OK\n");
	self->alnOpt.iHMMInputFiles = size;
	if(size ==0 )
		self->alnOpt.ppcHMMInput = NULL;
	else{
		self->alnOpt.ppcHMMInput = PyMem_Malloc(sizeof(char*) *  size);
        memset(self->alnOpt.ppcHMMInput, 0, sizeof(char*) *  size);
		for(i = 0; i < size; i++){
    			if(uti_PyString_to_pchar(&(self->alnOpt.ppcHMMInput[i]), PyTuple_GetItem(value,i)) == -1)
        			return -1;
		}
	}
return 0;
}

static PyObject * AlignOptions_getHMMInput(AlignOptions *self, void *closure){
    int i;
    PyObject *result;
    int size = self->alnOpt.iHMMInputFiles;
    result = PyTuple_New(size);
    if(size > 0){
        for(i = 0; i < size; i++){
		PyTuple_SetItem(result, i, uti_pchar_toPyString(self->alnOpt.ppcHMMInput[i]));
	}
    }
return result;
}







static PyGetSetDef AlignOptions_getseters[] ={
    {"distmatInfile", 
        (getter)AlignOptions_getDistmatInfile, 
        (setter)AlignOptions_setDistmatInfile,
        "distance matrix input file", NULL},
    {"distmatOutfile", 
        (getter)AlignOptions_getDistmatOutfile, 
        (setter)AlignOptions_setDistmatOutfile,
        "distance matrix output file", NULL},
    {"clustfile", 
        (getter)AlignOptions_getClustfile, 
        (setter)AlignOptions_setClustfile,
        "file with clustering information", NULL},
    {"guidetreeOutfile", 
        (getter)AlignOptions_getGuidetreeOutfile, 
        (setter)AlignOptions_setGuidetreeOutfile,
        "guidetree output file", NULL},
    {"guidetreeInfile", 
        (getter)AlignOptions_getGuidetreeInfile, 
        (setter)AlignOptions_setGuidetreeInfile,
        "guidetree input file", NULL},
    {"HMMInput", 
        (getter)AlignOptions_getHMMInput, 
        (setter)AlignOptions_setHMMInput,
        "HMM input files. index range: 0..iHMMInputFiles", NULL},
    {NULL}
};

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
    AlignOptions_methods,       /*tp_methods*/
    AlignOptions_members,       /*tp_members*/
    AlignOptions_getseters,     /*tp_getset*/
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
    {"Align", (PyCFunction) pyclustal_Align, METH_VARARGS,
    "Align(filename, alnOpt) -> result \n Same as getAlign, but with customer alnOpt(AlginOptions)"},
    {NULL, NULL, 0, NULL}
};

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
