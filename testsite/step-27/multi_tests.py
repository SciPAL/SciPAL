#!/usr/bin/env python2
import os
import sys
import string
from subprocess import call
import signal
import string
import time

STARTPATH = os.path.expanduser("~") + '/Praktikum_2013/testsite/step-28-build/'

Step28PrmPath = STARTPATH + 'prm/step-28.prm'
NumericsPrmPath = STARTPATH + 'prm/step-28_drs_numerics.prm'
TimingDir = STARTPATH + 'timing/'
#TimingDir = '/scratch/c.holme/timing/'
ExecPath = STARTPATH + 'step-28'

TimingPrefix = 'timing'


archList                = ['cuda', 'cpu']
n_threadList            = [12, 1]
BEM_Quad_OrderList      = [16, 12, 8, 4]
FE_Mapping_DegreeList   = [3, 2, 1]
FE_DegreeList           = [3, 2, 1]
N_init_refinementsList  = [1]
N_mesh_refinementsList  = [5]

###################################################################################################
def set_numerics_param(bem_quad_order, fe_mapping_degree, fe_degree,
        n_init_refinements, n_mesh_refinements, n_threads):
	prm = open(NumericsPrmPath)
	prm_list = prm.readlines()
	prm.close()
	i = 0
	for line in prm_list:	
		if(str.find(line,"set BEM Quad Order") != -1):
			prm_list[i] = "  set BEM Quad Order     = "+str(bem_quad_order)+"\n"
		if(str.find(line,"set FE Degree") != -1):
			prm_list[i] = "  set FE Degree  = "+str(fe_degree)+"\n"
		if(str.find(line,"set FE Mapping Degree") != -1):
			prm_list[i] = "  set FE Mapping Degree      = " + str(fe_mapping_degree)+"\n"
		if(str.find(line,"set N init refinements") != -1):
			prm_list[i] = "  set N init refinements = "+str(n_init_refinements)+"\n"
		if(str.find(line,"N mesh refinements") != -1):
			prm_list[i] = "  set N mesh refinements   = "+str(n_mesh_refinements)+"\n"
		if(str.find(line,"N threads") != -1):
			prm_list[i] = "  set N threads          = "+str(n_threads)+"\n"

		i +=1

	prm = open(NumericsPrmPath, "w")
	prm.writelines(prm_list)
	prm.close()

def set_arch(arch, bem_quad_order, fe_mapping_degree, 
                fe_degree, n_init_refinements, n_mesh_refinements):
# output_dir bauen
    outputDir = ''.join([ TimingPrefix, "_BEMQ-", str(bem_quad_order),"_FEMapp-", str(fe_mapping_degree)\
    ,"_FEDeg-", str(fe_degree)\
    ,"_InitRef-",str(n_init_refinements),"_MaxRef-",str(n_mesh_refinements)])


    prm = open(Step28PrmPath, "r")
    prm_list = prm.readlines()
    prm.close()
    i = 0
    for line in prm_list:
        if(str.find(line, "set Architecture") != -1):
            prm_list[i] = "   set Architecture        = " + str(arch) + "\n"
        if(str.find(line, "set Run directory") != -1):
            prm_list[i] = "     set Run directory    = " + str(TimingDir) + "\n"
        if(str.find(line, "set TimingOutputDirectory") != -1):
            prm_list[i] = "     set TimingOutputDirectory  = " \
                        + outputDir + "\n"

        i += 1

    prm = open(Step28PrmPath, "w")
    prm.writelines(prm_list)
    prm.close()

#for IMAGE in IMAGES:
#	for ALGORITHM in range(1,3):
#		for PRECISION in range(1,3):
#			set_param(IMAGE, ALGORITHM, PRECISION, ITERATIONS)
#			os.system("./step-16")

def main():
    for bem_quad_order in BEM_Quad_OrderList:      
        for fe_mapping_degree in FE_Mapping_DegreeList:
            for fe_degree in FE_DegreeList:
                for n_init_refinements in N_init_refinementsList:
                    for n_mesh_refinements in N_mesh_refinementsList:
                        for arch in archList:
                            for n_threads in n_threadList:
                                set_arch(arch, bem_quad_order, fe_mapping_degree, 
                                        fe_degree, n_init_refinements, n_mesh_refinements)
                                set_numerics_param(bem_quad_order, fe_mapping_degree, 
                                        fe_degree, n_init_refinements,
                                        n_mesh_refinements, n_threads)
#                        call(["/bin/cat", NumericsPrmPath ])
#                        call(["/bin/cat", Step28PrmPath ])
                                print "--------------------Timing " + arch +   " ----------------------"
                                os.system(ExecPath)
    plotAll()

################################## Plotting ###################################
PlotFileNames = ['step-28_cpu_timing.dat', 'step-28_cuda_timing.dat']

def plot_results(plotType, dir):
# header: 
#    1       2          3           4                  5         6          7                   8          9
# cycle n_dofs_per_be n_cells n_bem_quad_points n_bem_points FE_Degree FE_Mapping_Degree avg_iter_time total_Time

    file = open(dir + "/" + PlotFileNames[1])
    header = file.readline().split()
    params = file.readline().split()
    file.close
    numCols = len(header)# - 1
#    print numCols

    paramsDict = dict(zip(header, params))
    title =  str(plotType)
    title += ": #BEMQPoints: " + str(paramsDict['n_bem_quad_points'])
    title += ", #BEMPoints: " + str(paramsDict['n_bem_points'])
    title += ", #Cells: " + str(paramsDict['n_cells'])
    title += ", #DOFsPerBE: " + str(paramsDict['n_dofs_per_be'])


    plot =  "set term pdfcairo enhanced\n"
    plot += "set output \'" + dir + "/" + str(plotType) + ".pdf\'\n"
    plot += "set title \'" +  str(title) + "\'\n"

    plot += "set xlabel \'#cells\'\n"

    plot += "set logscale x\n"

    plot += "set ylabel \'Speedup time_{CPU}/time_{CUDA}\'\n"
    plot += "set ylabel \'Speedup time_{CPU}/time_{CUDA}\'\n"
    #plot += "set key top left\n"
    plot += "set key left Left reverse\n"
    plot += "file = \'< paste " + str(dir + "/" + PlotFileNames[0]) + " " + str(dir + "/" + PlotFileNames[1]) + " \'\n"
    plot += "numCols = " + str(numCols) + "\n"
    # With first numThreads
    plot += "plot file using 4:(column(10)/column(10+numCols)) every ::0::"
    plot += str(N_mesh_refinementsList[0]-1) +" w lp t \'total time with 1 thread\' pt 7 lc 1, \\\n"

    plot += "file using 4:(column(9)/column(9+numCols)) every ::0::"
    plot += str(N_mesh_refinementsList[0]-1) +" w lp t \'avg. iteration time with 1 thread\' lc 3 pt 13, \\\n"
    # With second numThreads
    plot += "file using 4:(column(10)/column(10+numCols)) every ::"+str(N_mesh_refinementsList[0]) +"::"
    plot += str(N_mesh_refinementsList[0]*2)+ " w lp t \'total time with 12 threads\' pt 7 lc 4, \\\n"

    plot += "file using 4:(column(9)/column(9+numCols)) every   ::"+str(N_mesh_refinementsList[0])+"::"
    plot += str(N_mesh_refinementsList[0]*2)+ " w lp t \'avg. iteration time with 12 threads\' lc 11 pt 13\n"



    plot += "reset"

#    print(''.join(["echo \"", plot, "\" | gnuplot "]))
    os.system(''.join(["echo \"", plot, "\" | gnuplot "]))




import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from re import compile,match


def plot_results_matplotlib(dir):
# header: 
# 0            1       2          3           4                  5         6          7                   8          9             
# cycle  n_threads n_dofs_per_be n_cells n_bem_quad_points n_bem_points FE_Degree FE_Mapping_Degree avg_iter_time total_Time
    # Read in header
    file = open(dir + "/" + PlotFileNames[1])
    header = file.readline().split()
    params = file.readline().split()
    file.close
    numCols = len(header)
    paramsDict = dict(zip(header, params))
    MT_threads=paramsDict['n_threads']

    #Read data from header
#    outputDir = ''.join([ TimingPrefix, "_BEMQ-", str(bem_quad_order),"_FEMapp-", str(fe_mapping_degree)\
#    ,"_FEDeg-", str(fe_degree)\
#    ,"_InitRef-",str(n_init_refinements),"_MaxRef-",str(n_mesh_refinements)])

    re_bemq = compile(r'.*BEMQ-(\d+)_.*')
    match_bemq = re_bemq.match(dir)
    if match_bemq:
        bem_quad_order = match_bemq.groups(1)[0]
    else:
        bem_quad_order = 0


    # Read in data
    cpu_mt = []
    cpu_st = []

    cuda_mt = []
    cuda_st = []

    file = open(dir + "/" + PlotFileNames[0])
    file.readline() #Jump over header
    for i in range(N_mesh_refinementsList[0]):
        cpu_mt.append(file.readline().split())
    file.readline() #jump over 2nd header
    for i in range(N_mesh_refinementsList[0]):
        cpu_st.append(file.readline().split())
    file.close()

    file = open(dir + "/" + PlotFileNames[1])
    file.readline() #Jump over header
    for i in range(N_mesh_refinementsList[0]):
        cuda_mt.append(file.readline().split())
    file.readline() #jump over 2nd header
    for i in range(N_mesh_refinementsList[0]):
        cuda_st.append(file.readline().split())
    file.close()


    #Construct the x-axis values, that means the n_cells
    xs = []
    for para_line in cuda_mt:
        xs.append(para_line[3])

    #Construct the y-value arrays, i.e. the speedups
    speedup_total_mt = []
    speedup_avg_mt = []
    speedup_total_st = []
    speedup_avg_st = []
    speedup_avg_mt_vs_st = []
    speedup_total_mt_vs_st = []

    for i in range(N_mesh_refinementsList[0]):
        speedup_total_mt.append(float(cpu_mt[i][9])/float(cuda_mt[i][9]))
        speedup_avg_mt.append(float(cpu_mt[i][8])/float(cuda_mt[i][8]))

        speedup_total_st.append(float(cpu_st[i][9])/float(cuda_st[i][9]))
        speedup_avg_st.append(float(cpu_st[i][8])/float(cuda_st[i][8]))

        speedup_total_mt_vs_st.append(float(cpu_st[i][9])/float(cuda_mt[i][9]))
        speedup_avg_mt_vs_st.append(float(cpu_st[i][8])/float(cuda_mt[i][8]))
        
        #Just for fun
#        speedup_total_mt_vs_st.append(float(cpu_mt[i][9])/float(cuda_st[i][9]))
#        speedup_avg_mt_vs_st.append(float(cpu_mt[i][8])/float(cuda_st[i][8]))


    fig = plt.figure(figsize=(10,6), dpi=200)
    ax = fig.add_subplot(111)
    
    ax.set_xlabel("Number of cells")
    ax.set_xscale("log")

    ax.set_ylabel("Speedup  "+ r'$t_{CPU}/t_{CUDA}$')

    LW=3.0

    titlestring ="Average iteration and total cell loop speedups" 
    titlestring += "\n #BEM quad points: " + str(paramsDict['n_bem_quad_points'])
    titlestring += ", #BEM points: " + str(paramsDict['n_bem_points'])
#    titlestring += ", #cells: " + str(paramsDict['n_cells'])
    titlestring += ", #DOFs per BE: " + str(paramsDict['n_dofs_per_be'])

    titlestring += "\n BEM quad order: " + str(bem_quad_order)
    titlestring += ", FE degree: " + str(paramsDict['FE_Degree'])
    titlestring += ", FE mapping degree: " + str(paramsDict['FE_Mapping_Degree'])
    plt.suptitle(titlestring)

    font = {'family' : 'monospace',
            'weight' : 'medium',
            'size'   : 11}
    matplotlib.rc('font', **font)
# Warning: ionly 12 threads allowd
    ax.plot(xs, speedup_total_mt, 'r-o', label="Total      12 threads", linewidth=LW)
    ax.plot(xs, speedup_avg_mt, 'r--d',  label="Avg. iter. 12 threads", linewidth=LW)

    ax.plot(xs, speedup_total_st, 'g-o', label="Total      1 thread", linewidth=LW)
    ax.plot(xs, speedup_avg_st, 'g--d',  label="Avg. iter. 1 thread", linewidth=LW)
    
    ax.plot(xs, speedup_total_mt_vs_st, 'b-o', label="Total     12 vs. 1 threads", linewidth=LW)
    ax.plot(xs, speedup_avg_mt_vs_st, 'b--d',  label="Avg. iter 12 vs. 1 threads", linewidth=LW)

    # Shink current axis's height by 15% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.15,
                             box.width, box.height * 0.80])

    plt.legend(loc='lower center', shadow=True, fancybox=True, ncol=3,
            bbox_to_anchor=(0.5, -0.33))

    ax.grid()


#    plt.show()
    outfile = dir + '/speedup_' + dir + '.pdf'
#    print outfile

    pp = PdfPages(outfile)
    pp.savefig(fig)
    pp.close()


PlotTypes = ['speedup']


def plotAll():
    os.chdir(TimingDir)

#    dir='timing_BEMQ-8_FEMapp-3_FEDeg-3_InitRef-1_MaxRef-3'
#    plot_results_matplotlib(dir)
#    sys.exit("blabla")


    for dir in [file for file in os.listdir(".") if (os.path.isdir(file) and
            file.startswith(TimingPrefix))]:
        for plotType in PlotTypes:
            print "Plotting in: " + dir
#            plot_results(plotType, dir)
            plot_results_matplotlib(dir)
            
        
def handler(signum, frame):
    sys.exit("Terminating...")



if __name__ == "__main__":
    signal.signal(signal.SIGTERM, handler)
    signal.signal(signal.SIGQUIT, handler)
    signal.signal(signal.SIGHUP , handler)
    signal.signal(signal.SIGABRT, handler)
    signal.signal(signal.SIGINT,  handler)
    if len(sys.argv) == 1:
        main()
    else:
        plotAll()

