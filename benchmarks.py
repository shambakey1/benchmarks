'''
Created on Sep 19, 2017

@author: Mohammed Elshambakey
'''

import os

def checkItemsInFile(fin,strin):
    ''' Check if a specific file has all specific items in a list '''

    with open(fin,"r") as fp:
        for line in fp:
            for str_item in strin:
                if str_item in line:
                    strin.remove(str_item)
                    if not strin:
                        return True    # All list items have been found in the file
    return False

def getFileNames(fin,sep):
    ''' Return a list of file names from a file containing a list of file paths '''

    out=[]
    with open(fin,"r") as f:
        for line in f:
            if sep in line:
                out.append(line.rsplit(sep,1)[1])
            else:
                out.append(line)
    return out

def getDiffFiles(f1in,f2in,fout,sep):
    ''' Get lines in f1 that do not exist in f2 '''

    cf1=getFileNames(f1in,sep)
    cf2=getFileNames(f2in,sep)
    diff=set(cf1).difference(set(cf2))
    print("cf1: "+str(len(cf1))+", cf2: "+str(len(cf2))+", diff: "+str(len(cf1)-len(cf2)))
    print("set cf1: "+str(len(set(cf1)))+", set cf2: "+str(len(set(cf2))))
    print("diff: "+str(len(diff))+", set diff: "+str(len(set(diff))))
    with open(fout,"w") as pfout:    
        for line in diff:
            pfout.write(line)

def getWrongResFiles(res_path):
    ''' Identify .res files (i.e., result files for system responsiveness) with wrong results '''

    resp_file_cond=["DESIRED STATE"]    # Correctness criteria for .res files
    wrongfiles=[]
    for f in os.listdir(res_path):
        if f.endswith(".res"):    # Check .res files. The correctness condition will be existence of 'DESIRED STATE' statement in f
            if not checkItemsInFile(os.path.join(res_path,f),list(resp_file_cond)): # Second parameter uses list() to pass a copy of the correctness criteria condition list, not the original list.
                wrongfiles.append(f)
    return wrongfiles

def getWrongZGESVOutFiles(res_path):
    ''' Identify .out files (i.e., result file for running zgesv lapacke benchmark code inside each container) with wrong results '''

    out_file_cond=["start_time","Entry Matrix","Right Rand Side","LAPACKE_zgesv","Solution","real","user","sys","end_time"]    # Correctness criteria for .res files
    wrongfiles=[]
    for f in os.listdir(res_path):
        if f.endswith(".out"):    # Check .res files. The correctness condition will be existence of 'DESIRED STATE' statement in f
            if not checkItemsInFile(os.path.join(res_path,f),list(out_file_cond)): # Second parameter uses list() to pass a copy of the correctness criteria condition list, not the original list.
                wrongfiles.append(f)
    return wrongfiles

def gen_zgesv_benchmarks(path,N_min,N_max,NRHS_min,NRHS_max):
    ''' Generate Benchmark files (c files) that use ZGESV function '''

    import random

    for N in range(N_min,N_max):
        for NRHS in range(NRHS_min,NRHS_max):
            fout=os.path.join(path,str(N)+"_"+str(NRHS)+".c")
            with open(fout,"w") as fp:
                fp.write("#include <stdlib.h>\n")
                fp.write("#include <stdio.h>\n")
                fp.write('#include "lapacke.h"\n\n')
                fp.write('extern void print_matrix( char* desc, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda);\n')
                fp.write('extern void print_int_vector( char* desc, lapack_int n, lapack_int* a);\n\n')
                fp.write('#define N '+str(N)+'\n')
                fp.write('#define NRHS '+str(NRHS)+'\n')
                fp.write('#define LDA N\n')
                fp.write('#define LDB NRHS\n\n')
                fp.write('int main() {\n')
                fp.write('\tlapack_int n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;\n')
                fp.write('\tlapack_int ipiv[N];\n')
                fp.write('\tlapack_complex_double a[LDA*N];\n')
                fp.write('\tlapack_complex_double b[LDB*N];\n\n')
                for i in range(0,N*N):
                    fp.write('\ta['+str(i)+']=lapack_make_complex_double('+str(random.uniform(-100,100))+','+str(random.uniform(-100,100))+');\n')
                fp.write('\n')
                for i in range(0,N*NRHS):
                    fp.write('\tb['+str(i)+']=lapack_make_complex_double('+str(random.uniform(-100,100))+','+str(random.uniform(-100,100))+');\n')
                fp.write("\n")
                fp.write('\tprint_matrix("Entry Matrix A", n, n, a, lda);\n')
                fp.write('\tprint_matrix("Right Rand Side", n, nrhs, b, ldb);\n')
                fp.write('\tprintf("\\n" );\n')
                fp.write('\tprintf( "LAPACKE_zgesv (row-major, high-level) Example Program Results\\n" );\n')
                fp.write('\tinfo = LAPACKE_zgesv( LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb );\n')
                fp.write('\tif( info > 0 ) {\n')
                fp.write('\t\tprintf("The diagonal element of the triangular factor of A,\\n" );\n')
                fp.write('\t\tprintf("U(%i,%i) is zero, so that A is singular;\\n", info, info );\n')
                fp.write('\t\tprintf("the solution could not be computed.\\n" );\n')
                fp.write('\t\texit( 1 );\n')
                fp.write('\t}\n')
                fp.write('\tprint_matrix( "Solution", n, nrhs, b, ldb);\n')
                fp.write('\tprint_matrix( "Details of LU factorization", n, n, a, lda);\n')
                fp.write('\tprint_int_vector( "Pivot indices", n, ipiv);\n')
                fp.write('\texit( 0 );\n')
                fp.write('}\n\n')
                fp.write('void print_matrix( char* desc, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda) {\n')
                fp.write('\tlapack_int i, j;\n')
                fp.write('\tprintf("\\n %s\\n", desc );\n')
                fp.write('\tfor( i = 0; i < m; i++ ) {\n')
                fp.write('\t\tfor( j = 0; j < n; j++ )\n')
                fp.write('\t\t\tprintf(" (%6.2f,%6.2f)", lapack_complex_double_real(a[i*lda+j]), lapack_complex_double_imag(a[i*lda+j]) );\n')
                fp.write('\t\tprintf("\\n" );\n')
                fp.write('\t}\n')
                fp.write('}\n\n')
                fp.write('void print_int_vector( char* desc, lapack_int n, lapack_int* a) {\n')
                fp.write('\tlapack_int j;\n')
                fp.write('\tprintf("\\n %s\\n", desc );\n')
                fp.write('\tfor( j = 0; j < n; j++ ) printf(" %6i", a[j] );\n')
                fp.write('\tprintf("\\n" );\n')
                fp.write('}\n')

def detRemZGESVTests(path,rows_min,rows_max,cols_min,cols_max,replicas_min,replicas_max,repeat_min,repeat_max):
    ''' Determine remaining lapacke zgesv tests '''

    import glob

    res_mis=[]
    out_mis=[]
    res_exist=glob.glob(os.path.join(path,"*.res"))
    out_exist=glob.glob(os.path.join(path,"*.out"))
    for r in range(rows_min,rows_max+1):
        for c in range(cols_min,cols_max+1):
            for repl in range(replicas_min,replicas_max+1):
                for rept in range(repeat_min,repeat_max+1):
                    out_f_pat=str(r)+"_"+str(c)+"_"+str(repl)+"_"+str(rept)    # Result file name pattern. Used to check existence of the result file
                    if (rept==repeat_min or rept==repeat_max) and not os.path.isfile(os.path.join(path,out_f_pat+".res")):
                        res_mis.append({"rows":r,"cols":c,"replicas":repl,"rept":rept})
                    if len(glob.glob(os.path.join(path,out_f_pat+"_*.out"))) != repl:
                        out_mis.append({"rows":r,"cols":c,"replicas":repl,"rept":rept})
    result_stat={"res_exist":res_exist,"res_mis":res_mis,"out_exist":out_exist,"out_mis":out_mis}
    return result_stat

def readInparam(fin):
    ''' Read required parameters from the input file "fin".
        fin: Input file with required parameters line by line. Each line is comma seperated
        Output is a complete list of parameters
    '''

    inparam=[]
    with open(fin,"r") as f:
        for line in f:
            inparam.append(line.strip().split(","))
    return inparam
        

def runBenchmarkTests(test=None,image_name="shambakey1/lapacke_bench",restart="on-failure",inparam=None):
    ''' Run benchmark tests according to the input list parameters (e.g., in addition to the common parameters like image name, and service name, additional benchmark parameters are provides like matrices specifications (rows, columns), number of service replicas, repeat number for each test ... etc). 
test: The required benchmark test (e.g., lapacke zgesv test). Each test can have different parameters
image_name: Image name at Docker hub with required libraries for the test. Default is shambakey1/lapacke_bench, but it can change for different benchmarks
'''

    import docker, sys

    if not test:
        print "Please enter a test"
        sys.exit()
    if not inparam:
        print "Not valid input parameters"
        sys.exit()
    client = docker.from_env()
    wrk_dir="/home"                # Work directory inside the running container. Currently, it is alawys HOME directory
    if test=="zgesv":            # If required test is lapacke zgesv
        mnts=[]                # List of mounts to be passed to created services
        wrk_dir_src="/home/ubuntu/benchmarks/zgesv"
        bench_com="./bench_task.sh"                # The command to be executed in each service container
        mnts.append(wrk_dir_src+":/home")            # Mount the zgesv directory inside service
        repeat_min=1
        repeat_max=10
        for inparam_item in inparam:
            env_list=[]                    # List of environment variables to be passed to created services
            rows=inparam_item[0]            # Matrix number of rows for current test
            cols=inparam_item[1]            # Matrix number of columns for current test
            replicas=inparam_item[2]        # Matrix number of replicas for current test
            if replicas>1:
                mode_type=docker.types.services.ServiceMode('replicated',replicas)
            else:
                mode_type=None
            rept=inparam_item[3]            # Matrix repeatition id for current test
            serv_name=str(rows)+"_"+str(cols)+"_"+str(replicas)+"_"+str(rept)    # Assign service name
            FIN=str(rows)+"_"+str(cols)            # The input C file to be executed in the container
            env_list.append("FIN="+FIN)
            FOUT="results/"+serv_name+"_{{.Task.ID}}"    # The outptut of the C benchmark file
            env_list.append("FOUT="+FOUT)
            client.services.create(image_name,bench_com,name=serv_name,workdir=wrk_dir,env=env_list,mounts=mnts,mode=mode_type,restart_policy=docker.types.services.RestartPolicy(condition=restart))
            if rept==repeat_min or rept==repeat_max:    # Check system responsiveness
                resp_cmd="time (docker service ps $(docker service ls -q)) &>> "+os.path.join(wrk_dir_src,"results",serv_name+".res")
                os.system(resp_cmd)

def zgesv_err_rem_files(res_path,rows_min,rows_max,cols_min,cols_max,replicas_min,replicas_max,repeat_min,repeat_max):
    ''' Remove wrong results and determine remaining experiments for the zgesv benchmark '''

    import glob
    
    err_pat_list=[]
    print "Determine wrong response files"    
    err_res_list=getWrongResFiles(res_path)                    # Determine wrong response files
    for f in err_res_list:
        err_pat_list.append(os.path.basename(f).rsplit("_",1)[0])    # Recored the pattern of wrong response files
    print "Determine wrong out files"
    err_out_list=getWrongZGESVOutFiles(res_path)                # Determine wrong out files
    for f in err_out_list:
        err_pat_list.append(os.path.basename(f).rsplit("_",2)[0])    # Record the pattern of wrong out files
    err_set=set(err_pat_list)
    print "Remove wrong response and out files"
    for f in err_set:
        fin=os.path.join(res_path,f+"_*")
        out=glob.glob(fin)
        for ferr in out:
            os.remove(ferr)            # Remove wrong responase and out file
    print "Determin remaining zgesv tests"
    return detRemZGESVTests(res_path,rows_min,rows_max,cols_min,cols_max,replicas_min,replicas_max,repeat_min,repeat_max)

def tempRemZGESVExp():
    ''' This is just temporary function during development. Do not use it '''

    res_path="/f/personal/postdoc/uncc/vifi/benchmarks/zgesv/results"
    rows_min=1
    rows_max=2
    cols_min=1
    cols_max=10
    replicas_min=1
    replicas_max=50
    repeat_min=1
    repeat_max=10
    out=zgesv_err_rem_files(res_path,rows_min,rows_max,cols_min,cols_max,replicas_min,replicas_max,repeat_min,repeat_max)
    out_final=set()
    for i in out['res_mis']:
        out_final.add((i['rows'],i['cols'],i['replicas']))
    for i in out['out_mis']:
        if (i['rows'],i['cols'],i['replicas']) not in out_final:
            out_final.add((i['rows'],i['cols'],i['replicas']))
    with open("test_mis.out","w") as f:
        for i in out_final:
            for j in range(10):
                f.write(str(i[0])+","+str(i[1])+","+str(i[2])+","+str(j+1)+"\n")

def collectZGESVOutFileResults(fin):
    ''' Extract results from file with results of lapacke zgesv benchmark '''
   
inparam=readInparam("/home/ubuntu/benchmarks/mis_test_rows_1.out")
runBenchmarkTests(test='zgesv',image_name="shambakey1/lapacke_bench",restart="on-failure",inparam=inparam)
    

