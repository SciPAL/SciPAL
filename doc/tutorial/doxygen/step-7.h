/**
 * @page step_7                 The step-7 tutorial program
@htmlonly
<table class="tutorial" width="100%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
      </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#ClassWaveForms">Class: WaveForms</a>
      <ul>
        <li><a href="#Functiongenerate_waveforms">Function: generate_waveforms</a>
        <li><a href="#Functionnormalize_waveforms">Function: normalize_waveforms</a>
        <li><a href="#Functionoriginal_waveforms">Function: original_waveforms</a>
      </ul>
        <li><a href="#Declarationsofthekernelwrapperfunctions">Declarations of the kernel-wrapper-functions</a>
        <li><a href="#Kernel__subtract_array_from_matrix">Kernel: __subtract_array_from_matrix</a>
        <li><a href="#Kernel__MM_scalar">Kernel: __MM_scalar</a>
        <li><a href="#Kernel__Det">Kernel: __Det</a>
        <li><a href="#Kernel__P_ij">Kernel: __P_ij</a>
        <li><a href="#Kernel__u_ij">Kernel: __u_ij</a>
        <li><a href="#Kernel__z_ij">Kernel: __z_ij</a>
      <ul>
        <li><a href="#Functionreduction">Function: reduction</a>
      </ul>
        <li><a href="#Kernel__sum_array">Kernel: __sum_array</a>
        <li><a href="#Kernel__sum_log_array">Kernel: __sum_log_array</a>
        <li><a href="#Kernel__sum_array2D">Kernel: __sum_array2D</a>
        <li><a href="#Kernel__y">Kernel: __y</a>
        <li><a href="#Kernel__MxM">Kernel: __MxM</a>
        <li><a href="#Kernel__SMX">Kernel: __SMX</a>
      <ul>
        <li><a href="#Functionsubtract_array_from_matrix">Function:subtract_array_from_matrix</a>
        <li><a href="#FunctionMM_scalar">Function: MM_scalar</a>
        <li><a href="#FunctionDet">Function: Det</a>
        <li><a href="#Functionassemble_P_ij">Function: assemble_P_ij</a>
        <li><a href="#Functionassemble_u_ij">Function: assemble_u_ij</a>
        <li><a href="#Functionassemble_z_ij">Function: assemble_z_ij</a>
        <li><a href="#Functioncalculate_y">Function: calculate_y</a>
        <li><a href="#Functionsum_array">Function: sum_array</a>
        <li><a href="#Functionsum_matrix_cols">Function: sum_matrix_cols</a>
        <li><a href="#Functionsum_array2D">Function: sum_array2D</a>
        <li><a href="#FunctionMxM">Function: MxM</a>
        <li><a href="#FunctionSMX">Function: SMX</a>
      </ul>
        <li><a href="#ClassShohamEM">Class: ShohamEM</a>
        <li><a href="#ClassKernelTest">Class: KernelTest</a>
        <li><a href="#StructRightMSolve">Struct: RightMSolve</a>
      <ul>
        <li><a href="#Operator">Operator: = </a>
      </ul>
        <li><a href="#StructRightMSolveTr">Struct: RightMSolveTr</a>
      <ul>
        <li><a href="#Operator">Operator: = </a>
        <li><a href="#Operator">Operator: *</a>
        <li><a href="#Operator">Operator: *</a>
        <li><a href="#ConstructorShohamEM">Constructor: ShohamEM</a>
        <li><a href="#Functiondump">Function: dump</a>
        <li><a href="#Functionrun_Shoham">Function: run_Shoham</a>
        <li><a href="#Functioninitialize_Shoham">Function: initialize_Shoham</a>
        <li><a href="#Functione_step_Shoham">Function: e_step_Shoham</a>
        <li><a href="#Functionm_step_Shoham">Function: m_step_Shoham</a>
        <li><a href="#Functionpurge_step_Shoham">Function: purge_step_Shoham</a>
        <li><a href="#DestructorShohamEM">Destructor : ~ShohamEM</a>
        <li><a href="#ConstructorKernelTest">Constructor: KernelTest</a>
        <li><a href="#Functionrun">Function: run</a>
        <li><a href="#DestructorKernelTest">Destructor : ~KernelTest</a>
      </ul>
        <li><a href="#ClassSimParams">Class: SimParams</a>
        <li><a href="#ClassPCA">Class: PCA</a>
        <li><a href="#ClassMyFancySimulation">Class: MyFancySimulation</a>
      <ul>
        <li><a href="#ConstructorPCA">Constructor: PCA</a>
        <li><a href="#Functionrun">Function: run</a>
        <li><a href="#Functionfactorize">Function: factorize</a>
        <li><a href="#Functiondeclare">Function: declare</a>
        <li><a href="#Functionget">Function: get</a>
        <li><a href="#Functiongenerate_waveforms">Function: generate_waveforms</a>
        <li><a href="#ConstructorMyFancySimulation">Constructor: MyFancySimulation</a>
        <li><a href="#Funktionrun">Funktion: run</a>
        <li><a href="#Funktionrun_cluster_analysis">Funktion: run_cluster_analysis</a>
        <li><a href="#Funktionrun_simulation">Funktion: run_simulation</a>
        <li><a href="#Funktionrun_kerneltest">Funktion: run_kerneltest</a>
        <li><a href="#Funktionmain">Funktion: main</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
      </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
    <ul>
        <li><a href="#plain-ClassWaveForms">Class: WaveForms</a>
      <ul>
        <li><a href="#plain-Functiongenerate_waveforms">Function: generate_waveforms</a>
        <li><a href="#plain-Functionnormalize_waveforms">Function: normalize_waveforms</a>
        <li><a href="#plain-Functionoriginal_waveforms">Function: original_waveforms</a>
      </ul>
        <li><a href="#plain-Declarationsofthekernelwrapperfunctions">Declarations of the kernel-wrapper-functions</a>
        <li><a href="#plain-Kernel__subtract_array_from_matrix">Kernel: __subtract_array_from_matrix</a>
        <li><a href="#plain-Kernel__MM_scalar">Kernel: __MM_scalar</a>
        <li><a href="#plain-Kernel__Det">Kernel: __Det</a>
        <li><a href="#plain-Kernel__P_ij">Kernel: __P_ij</a>
        <li><a href="#plain-Kernel__u_ij">Kernel: __u_ij</a>
        <li><a href="#plain-Kernel__z_ij">Kernel: __z_ij</a>
      <ul>
        <li><a href="#plain-Functionreduction">Function: reduction</a>
      </ul>
        <li><a href="#plain-Kernel__sum_array">Kernel: __sum_array</a>
        <li><a href="#plain-Kernel__sum_log_array">Kernel: __sum_log_array</a>
        <li><a href="#plain-Kernel__sum_array2D">Kernel: __sum_array2D</a>
        <li><a href="#plain-Kernel__y">Kernel: __y</a>
        <li><a href="#plain-Kernel__MxM">Kernel: __MxM</a>
        <li><a href="#plain-Kernel__SMX">Kernel: __SMX</a>
      <ul>
        <li><a href="#plain-Functionsubtract_array_from_matrix">Function:subtract_array_from_matrix</a>
        <li><a href="#plain-FunctionMM_scalar">Function: MM_scalar</a>
        <li><a href="#plain-FunctionDet">Function: Det</a>
        <li><a href="#plain-Functionassemble_P_ij">Function: assemble_P_ij</a>
        <li><a href="#plain-Functionassemble_u_ij">Function: assemble_u_ij</a>
        <li><a href="#plain-Functionassemble_z_ij">Function: assemble_z_ij</a>
        <li><a href="#plain-Functioncalculate_y">Function: calculate_y</a>
        <li><a href="#plain-Functionsum_array">Function: sum_array</a>
        <li><a href="#plain-Functionsum_matrix_cols">Function: sum_matrix_cols</a>
        <li><a href="#plain-Functionsum_array2D">Function: sum_array2D</a>
        <li><a href="#plain-FunctionMxM">Function: MxM</a>
        <li><a href="#plain-FunctionSMX">Function: SMX</a>
      </ul>
        <li><a href="#plain-ClassShohamEM">Class: ShohamEM</a>
        <li><a href="#plain-ClassKernelTest">Class: KernelTest</a>
        <li><a href="#plain-StructRightMSolve">Struct: RightMSolve</a>
      <ul>
        <li><a href="#plain-Operator">Operator: = </a>
      </ul>
        <li><a href="#plain-StructRightMSolveTr">Struct: RightMSolveTr</a>
      <ul>
        <li><a href="#plain-Operator">Operator: = </a>
        <li><a href="#plain-Operator">Operator: *</a>
        <li><a href="#plain-Operator">Operator: *</a>
        <li><a href="#plain-ConstructorShohamEM">Constructor: ShohamEM</a>
        <li><a href="#plain-Functiondump">Function: dump</a>
        <li><a href="#plain-Functionrun_Shoham">Function: run_Shoham</a>
        <li><a href="#plain-Functioninitialize_Shoham">Function: initialize_Shoham</a>
        <li><a href="#plain-Functione_step_Shoham">Function: e_step_Shoham</a>
        <li><a href="#plain-Functionm_step_Shoham">Function: m_step_Shoham</a>
        <li><a href="#plain-Functionpurge_step_Shoham">Function: purge_step_Shoham</a>
        <li><a href="#plain-DestructorShohamEM">Destructor : ~ShohamEM</a>
        <li><a href="#plain-ConstructorKernelTest">Constructor: KernelTest</a>
        <li><a href="#plain-Functionrun">Function: run</a>
        <li><a href="#plain-DestructorKernelTest">Destructor : ~KernelTest</a>
      </ul>
        <li><a href="#plain-ClassSimParams">Class: SimParams</a>
        <li><a href="#plain-ClassPCA">Class: PCA</a>
        <li><a href="#plain-ClassMyFancySimulation">Class: MyFancySimulation</a>
      <ul>
        <li><a href="#plain-ConstructorPCA">Constructor: PCA</a>
        <li><a href="#plain-Functionrun">Function: run</a>
        <li><a href="#plain-Functionfactorize">Function: factorize</a>
        <li><a href="#plain-Functiondeclare">Function: declare</a>
        <li><a href="#plain-Functionget">Function: get</a>
        <li><a href="#plain-Functiongenerate_waveforms">Function: generate_waveforms</a>
        <li><a href="#plain-ConstructorMyFancySimulation">Constructor: MyFancySimulation</a>
        <li><a href="#plain-Funktionrun">Funktion: run</a>
        <li><a href="#plain-Funktionrun_cluster_analysis">Funktion: run_cluster_analysis</a>
        <li><a href="#plain-Funktionrun_simulation">Funktion: run_simulation</a>
        <li><a href="#plain-Funktionrun_kerneltest">Funktion: run_kerneltest</a>
        <li><a href="#plain-Funktionmain">Funktion: main</a>
      </ul>
      </ul>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Introduction"></a><h1>Introduction</h1>

<p> 
To automatically classify the different shapes in the clustering problem, which we are facing in this project,
we have the option to use the Gaussian mixture model algorithm to determine the cluster responsibilities.
Although "recent studies have shown that the Gaussian model does not accurately capture the multivariate statistic of the waveform samples' distribution" <a href="http://dx.doi.org/10.1016%2fS0165-0270(03)00120-1">(Shy Shoham et al. 2003)</a>.<br>
In his publication "Robust, automatic spike sorting using mixtures of multivariate t-distributions" from 2003,
Shoham presents "further data demonstrating non-Gaussian statistic, and show that the multivariate t-distribution,
a wide-tailed family of distributions, provides a significantly better fit to the true statistic".<br>
Furthermore Shoham provides an algorithm which "is statistically plausible, simple and well-behaved and can effectively deal with many real data sets".
This algorithm is based on the mixture decomposition algorithm for multivariate t-distributions <a href="http://www.maths.uq.edu.au/%7Egjm/pm_sc00.pdf">(Peel and McLachlan, 2000)</a>
which requires computation of twice as many hidden variables as in Gaussian mixture algorithms and also involves
additional computation step for the 'degree of freedom' parameter.<br>
Instead of applying the algorithm directly, it would be applied in conjunction with an efficient model selection scheme developed by <a href="http://dx.doi.org/10.1109%2f34.990138">(Figueiredo and Jain, 2002)</a>,
which maximize a penalized log-likelihood with penalty based on the minimum message length criterion <a href="http://links.jstor.org/sici?sici=0035-9246%281987%2949%3A3%3C240%3AEAIBCC%3E2.0.CO%3B2-M">(Wallace and Freeman, 1987)</a>.
With this approach we get $L_{p}$ noted below, where $N$ is the number of parameters specifying each mixture component.[2]<br>

The complete algorithm goes as follows:<br>

<b>Algorithm:<br>
</p>
<p>
legend:<br>
</b>
</p>
<p>
$p :=$ space dimension.<br>
$g :=$ Number of cluster centers.<br>
$g_{min} :=$ Minimal number of cluster centers.<br>
$g_{max} :=$ Maximal number of cluster centers we are searching for.<br>
$n :=$ Number of data points.<br>
$N :=$ Number of parameters per mixture component.<br>
$\pi_j :=$ relative frequency of cluster $j$.<br>
$\Sigma_j :=$ Covariance matrix of the cluster $j$.<br>
$\Delta_j := \sqrt{det(\Sigma_j)}$<br>
$v :=$ Degree of freedom (DOF) parameter.<br>
$L :=$ log-likelihood function.<br>

</p>
<p>
<b>Initialization:</b></br>
<ul>
<li>use clustering method i.e K-means to determine centers $\vec{\mu}_0 ,..., \vec{\mu}_{g-1}$</li>
  <li>$g = g_{max}$</li>
  <li>Set</li>
  <ul>
    <li>$\pi_1 ,..., \pi_g = \frac{1}{g}$, $\pi \in \mathbb{K}^{g}$</li>
    <li>$\Sigma_1 ,..., \Sigma_g = I$, $\Sigma_j \in \mathbb{R}^{g \times g}$, $\Sigma \in \left (\mathbb{R}^{g \times g} \right )^{g}$</li>
    <li>$v = 50$</li>
    <li>$P_{ij}\in \mathbb{R}^{nxg}$</li>
    <li>$L_{max} = -\infty$</li>
    <li>$N$</li>
    <li>$x = \left \{ \vec{x}_0,...,\vec{x}_{n-1} \right \}$, $x \in \mathbb{R}^{p \times n}$</li>
    <li>$\mu = \left \{ \vec{\mu}_0,...,\vec{\mu}_{g-1} \right \}$, $\mu \in \mathbb{R}^{p \times g}$</li>
    <li>
    $\delta_{ij}:=\delta(\vec{x}_{i},\vec{\mu}_{j};\Sigma_{j})=(\vec{x}_{i}-\vec{\mu}_{j})^{T} I^{-1}(\vec{x}_{i}-\vec{\mu}_{j})$
    </li>
  </ul>
</ul>
</p><br />
<b>Iteration:</b></br>
<p>
<b>while</b> $g \geq g_{min}$<br><br />
&nbsp; &nbsp; <b>E-Step</b>
<ul>
<li>weighted squared distance between i-th datapoint and j-th cluster mean</li>
<li>(Mahalanobis distance)</li>
</ul>
\f{eqnarray}\delta_{ij}:=\delta(\vec{x}_{i},\vec{\mu}_{j};\Sigma_{j})&=&(\vec{x}_{i}-\vec{\mu}_{j})^{T}\Sigma_{j}^{-1}(\vec{x}_{i}-\vec{\mu}_{j})
    \f}<br />
<ul>
<li>likelihood of i-th datapoint assigned to j-th cluster</li>
</ul>
    \f{eqnarray}P_{ij}&=&\frac{\Gamma(\frac{v+p}{2})}{\Gamma(\frac{v}{2})(\pi v)^{p/2}\Delta_j}\frac{1}{(1+\frac{\delta_{ij}}{v})^{(v+p)/2}}
    \f}<br>
<ul>
<li>update $L$ with</li>
</ul>
    \f{eqnarray}L_p = \sum_{i=1}^n \log \sum_{j=1}^g P_{ij}\pi_j -\left [ \frac{N}{2}\sum_{j=1}^{g} \log \frac{n\pi_j}{12}+\frac{g}{2}+\frac{g(N+1)}{2} \right ]
    \f}<br>
&nbsp; &nbsp; <b>if</b>$(\Delta L < 0.1 \& \Delta v < 10^{-2})$: convergence reached<br>
&nbsp; &nbsp; &nbsp; &nbsp; <b>return;</b><br />
<ul>
<li>update memberships</li>
</ul>
    \f{eqnarray}\hat{z}_{ij}&=&\frac{P_{ij} \pi_j}{\sum_{l=1}^{g} P_{il} \pi_l}
    \f}<br>
<ul>
<li>and the weights</li>
</ul>
    \f{eqnarray}\hat{u}_{ij}&=&\frac{p+v}{\delta_{ij}+v}
    \f}<br>
<ul>
<li>calculate $y$ wich is needed to determine the parameter $v$</li>
</ul>
<br>
    \f{eqnarray}
    y \equiv -\sum_{i=1}^{n}\sum_{j=1}^{g}\hat{z}_{ij}\left [\psi\left(\frac{p+v_{old}}{2} \right)+ \log \left(\frac{2}{\delta_{ij}+v_{old}} \right)-\hat{u}_{ij} \right ]/n
    \f}<br>

&nbsp; &nbsp; <b>M-Step</b><br>
&nbsp; &nbsp; &nbsp; &nbsp; <b>while</b>$|\sum_{j=1}^{g}\pi_j -1|>{10}^{-4}$<br>
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; <b>For</b>$j=1:g$<br>
<ul>
<li>Update $\pi_j$ with</li>
</ul>
    \f{eqnarray}
    \pi_{j}^{(k)}&=&\frac{max\left(\sum_{i=1}^{n}\frac{P_{ij}\pi_j^{(k-1)}}{\sum_{l=1}^{g}P_{il}\pi_l^{(k-1)}}-\frac{N}{2},0 \right )}{n-\frac{gN}{2}}
    \f}<br>
<ul>
<li>update $g$ with number of $\pi_j>0$</li>
</ul>
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; <b>End For</b><br>
&nbsp; &nbsp; &nbsp; &nbsp; <b>End While</b><br>
<ul>
    <li>Purge components where $\pi_j = 0$</li>
    <li>calculate new cluster means $\mu$ and their respective standard deviation $\Sigma$ with</li>
</ul>
    \f{eqnarray}\mu_j &=& \frac{\sum_{i=1}^{n}\hat{z}_{ij}\hat{u}_{ij}\vec{x}_{i}}{\sum_{i=1}^{n}\hat{z}_{ij}\hat{u}_{ij}}
    \f}<br>
    \f{eqnarray}\Sigma_j &=& \frac{\sum_{i=1}^{n}(\hat{z}_{ij}\hat{u}_{ij})(\vec{x}_i - \vec{\mu}_j)(\vec{x}_i - \vec{\mu}_j)^T}{\sum_{i=1}^{n}\hat{z}_{ij}\hat{u}_{ij}}
    \f}<br>
    <br>
<ul>
<li>using Cholesky decomposition leads to</li>
</ul>
\f{eqnarray}\Sigma_j=L^TL\f}<br />
    \f{eqnarray}\Delta_j &=& \sqrt{det(\Sigma_j)}=\prod_i L_{ii}
    \f}<br>
<ul>
<li>update $v$</li>
</ul>
    \f{eqnarray}v_{new}&=&\frac{2}{y+\log y-1}+0.0416 \left(1+erf \left(0.6594*\log \left(\frac{2.1971}{y+\log y-1} \right) \right) \right)
    \nonumber \\
    \f}<br>
&nbsp; &nbsp; &nbsp; &nbsp; <b>If</b> $L > L_{max}$<br>
<ul>
    <li>$L_{max}=L$</li>
    <li>Store parameters {$\pi,\mu,\Sigma$} as 'optimal'</li>
    <li>Set smallest component to zero (regarding to their $\pi$ value)</li>
    <li>$g=g-1$</li>
</ul>
&nbsp; &nbsp; &nbsp; &nbsp; <b>Else</b><br>
<ul>
    <li><b>break;</b></li>
</ul>
&nbsp; &nbsp; &nbsp; &nbsp; <b>End if</b><br>
<b>End While</b><br>
</p>

<p>
Literature:<br>
1. <a href="http://www.bm.technion.ac.il/Labs/niel/Public%20Data/Publications/Shoham_etal_JNM2003.pdf">Robust, automatic spike sorting using mixtures of multivariate t-distributions (Shy Shoham et al., 2003)</a><br>
2. <a href="http://www.lx.it.pt/~mtf/IEEE_TPAMI_2002.pdf ">Unsupervised Learning of Finite Mixture Models (Figueirdo and Jain, 2002)</a><br>
</p>

<p>
<b>Parallelzation structure:</b><br>

In this paragraph we discuss how we intend to implement each Kernel needed for the individual equations.<br>
<br>
<b><i>$P_{ij}$ Kernel equation</i></b><br>
first we split the equation in three parts:<br>
\f{eqnarray}P_{ij}&=&\frac{\Gamma(\frac{v+p}{2})}{\Gamma(\frac{v}{2})(\pi v)^{p/2}}\cdot\frac{1}{\Delta_j}\cdot\frac{1}{(1+\frac{\delta_{ij}}{v})^{(v+p)/2}}
\f}<br>
The first part: $\gamma := \frac{\Gamma(\frac{v+p}{2})}{\Gamma(\frac{v}{2})(\pi v)^{p/2}}$ should be pre-calculated using the CPU and would be handed out
to the Kernel as an argument.<br>
<br>
The second part: $\frac{1}{\Delta_j}$ with $\Delta_j$ being a pointer on a small array with $\sqrt{det \Sigma_j}$<br>
<br>
For the last part: $\frac{1}{(1+\frac{\delta_{ij}}{v})^{(v+p)/2}}$ each $\delta_{ij}$ would be loaded exactly once.<br>
<br>
Therefor our function parameter for the kernel are:<br>
<ul>
    <li> $\{\gamma, v , p\} \in \mathbb{R}$</li>
    <li> $n, g$ as space dimension</li>
    <li> const T * Delta $\in \mathbb{R}^{g}$</li>
    <li> const T * delta $\in \mathbb{R}^{n \times g}$</li>
    <li> T * const P_ij $\in \mathbb{R}^{n \times g}$</li>
</ul>

For the thread and Grid dimensions we choose<br>
<ul>
    <li> 1-dimensional threadblocks to handel the column sections</li>
    <li> 2-dimensional thread-Grids with:
        <ul>
            <li>gridDim.x = column index --> which $\Delta_j$</li>
            <li>gridDim.y = top row of the threadblock</li>
        </ul>
The indices $i$ and $j$ would be:<br>
    <li> $I = $ gridDim.y * blockDim.x + threadIdx.x</li>
    <li> $j = $ gridDim.x</li>
</ul>
Finally we calculate $P_{ij}$ and store it in the P_ij Matrix.<br>
</p>

<b><i>$\frac {\sum_{i=1}^{n}\hat{z}_{ij}\hat{u}_{ij}\vec{x}_i}{\sum_{i=1}^{n}\hat{z}_{ij}\hat{u}_{ij}}$ Kernel function</i></b><br>
first we split the equation in two parts:<br>
The denominator: $\sum_{i=1}^{n}\hat{z}_{ij}\hat{u}_{ij}$ is calculated with the sum_array2D Kernel.<br>
<br>
The numerator: $\sum_{i=1}^{n}\hat{z}_{ij}\hat{u}_{ij}\vec{x}_i$ present a challenge in loading the vector $x_i$ only once
and calculating the sum of the matrices which are created by the vector product of $\hat{z}_{ij}\hat{u}_{ij} \vec{x}_i$<br>
<br>
So we calculate the following:<br>
$\mu = \begin{bmatrix}
. & . & . & . & . \\
. & . & . & . & . \\
. & . & . & . & . \\
. & . & . & . & . \\
. & . & . & . & . \\
. & . & . & . & . \\
. & . & . & . & . \\
. & . & . & . & . \\
. & . & . & . & . \\
. & . & . & . & . \\
\end{bmatrix} = \begin{pmatrix}
.\\
.\\
.\\
\vec{x}_i\\
\in\\
\mathbb{R}^p\\
.\\
.\\
.\\

\end{pmatrix} \cdot (\hat{z}_{i0}\hat{u}_{i0},...,\hat{z}_{ig}\hat{u}_{ig})\in \mathbb{R}^g  + \begin{pmatrix}
.\\
.\\
.\\
\vec{x}_{i+1}\\
\in\\
\mathbb{R}^p\\
.\\
.\\
.\\

\end{pmatrix} \cdot (\hat{z}_{i+1 0}\hat{u}_{i+1 0},...,\hat{z}_{i+1g}\hat{u}_{i+1g})\in \mathbb{R}^g + ...$<br>

Therefore our function parameter for the kernel are:<br>
<ul>
    <li> $n, g, p$ as space dimension $\in \mathbb{Z}$</li>
    <li> const T * z_ij $\in \mathbb{R}^{n \times g}$</li>
    <li> const T * u_ij $\in \mathbb{R}^{n \times g}$</li>
    <li> const T * x $\in \mathbb{R}^{n \times p}$</li>
    <li> T * $\mu \in \mathbb{R}^{g \times p}$</li>
</ul>

For the threadblock and Grid dimensions we choose<br>
<ul>
    <li> 3-dimensional threadblocks with:</li>
        <ul>
            <li>blockDim.x = handle the $i$ components of the $x$ array.<br></li>
            <li>blockDim.y = the $g$ components of the $\hat{z}_{ij}\hat{u}_{ij}$<br></li>
            <li>blockDim.z = the $p$ components.<br></li>
            <li>since we want to start 512 threads the blocksize will be dim3(32,4,4);</li>
            <li>32 refers to the number of threads in one warp which will enable us to easily calculate the sum of the individual portions.</li>
            <li> since we have no option to get the space dimensions $p , g$ during the compilation time, we portion the number of data points in 4 data points in each direction</li>
        </ul>
    <li> 3-dimensional thread-Grids with:
        <ul>
            <li>gridDim.x = 32 </li>
            <li>gridDim.y = (g+4-1)/4 </li>
            <li>gridDim.z = (p+4-1)/4 </li>
        </ul>
</ul>
</p>
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * @code
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef WAVEFORMS_GENERATOR_H
 * #define WAVEFORMS_GENERATOR_H
 * 
 * #include <vector>
 * #include <time.h>
 * #include <iostream>
 * #include <vector>
 * #include <fstream>
 * #include <QDir>
 * #include <QString>
 * #include <QVector>
 * 
 * @endcode
 * 
 * boost random generator
 * 
 * @code
 * #include <boost/random/mersenne_twister.hpp>
 * #include <boost/random/normal_distribution.hpp>
 * #include <boost/random/variate_generator.hpp>
 * #include <boost/random/uniform_real.hpp>
 * #include <boost/math/distributions/students_t.hpp>
 * 
 * @endcode
 * 
 * 
 * <a name="ClassWaveForms"></a> 
 * <h3>Class: WaveForms</h3>
 * 
 * 
 * Generates the initial Test data
 * 
 * @code
 * class WaveForms{
 * 
 * public:
 *     WaveForms();
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functiongenerate_waveforms"></a> 
 * <h4>Function: generate_waveforms</h4>
 * This function generates the initial waveforms for testing the ShohamEM algorithm
 * 
 * 
 * @param row : number of rows (number of waves)
 * @param col : number of columns (space dimension)
 * @param waveforms : the matrix where the data are written
 * @param n_waves : number of waves from each form
 * @param S : noise intensity
 * @param waves_parameter : conataints the parameter needed for generating the diffrent waveforms
 * 
 * @code
 *     template <typename D>
 *     static void generate_waveforms(unsigned int row,
 *                                    unsigned int col,
 *                                    dealii::FullMatrix<D> &waveforms,
 *                                    const std::vector<int> &n_waves,
 *                                    D S,
 *                                    const std::vector<double> &waves_parameter)
 *     {
 * @endcode
 * 
 * random number generator, generates numbers between -1. and 1.
 * 
 * @code
 *         boost::mt19937 rng;
 *         boost::uniform_real<D> u(-1., 1.);
 *         boost::variate_generator<boost::mt19937&, boost::uniform_real<D> > gen(rng, u);
 * 
 * @endcode
 * 
 * random number generator, generates numbers between 0. and 1.
 * 
 * @code
 *         boost::mt19937 rng2;
 *         boost::uniform_real<D> v(0., 1.);
 *         boost::variate_generator<boost::mt19937&, boost::uniform_real<D> > gen2(rng2, v);
 * 
 * 
 *         double form;
 *         int form_id;
 *         double X, X1, X2, sig, sig1, sig2, nu, t, Y, noise, l, t0;
 *         double P0,P1,P2,P3,P4,P5,P6;
 * 
 * @endcode
 * 
 * store the waveforms parameter in  different variables
 * 
 * @code
 *         P0 = waves_parameter[0];
 *         P1 = waves_parameter[1];
 *         P2 = waves_parameter[2];
 *         P3 = waves_parameter[3];
 *         P4 = waves_parameter[4];
 *         P5 = waves_parameter[5];
 *         P6 = waves_parameter[6];
 * 
 * @endcode
 * 
 * for each row/datapoint determine which form to write by calculating the commonness of each form
 * 
 * @code
 *         for(int j = 0; j < row ; j++){
 * 
 *             form = gen2();
 * 
 *             if(form<(n_waves[0]/(double)row))
 *                 form_id=1;
 *             else{
 *                 if(form > ((n_waves[0]+n_waves[1])/(double)row))
 *                     form_id = 3;
 *                 else
 *                     form_id= 2;
 *             }
 * 
 * @endcode
 * 
 * using the switch instruction and the calculated form ID, decide which form is to be written in the j_th row.
 * 
 * @code
 *             switch(form_id){
 * @endcode
 * 
 * \f{equation} Y = \frac{1}{\tau} \ast \left(x-x_{0}\right)\exp-\left(\frac{x-x_{0}}{\sigma_{0}\sqrt{2}}\right)^2 + R \f}
 * 
 * @code
 *             case 1:{
 *                 noise = 1+(gen()*S);
 * 
 *                 X = P0*(noise);
 *                 sig = P1*(noise);
 *                 t = 0.1*(noise);
 *                 l = 0.;
 * 
 *                 for(int k = 0; k<col; k++){
 *                     Y = gen()*S;
 *                     l = k/double(col);
 *                     waveforms(j,k) = (1+Y)+(1/t)*((l-X)*exp(-(std::pow(((l-X)/(sig*sqrt(2))),2))));
 *                 }
 *                 break;
 *             }
 * @endcode
 * 
 * \f{equation} Y = cos^2(\frac{t-t_0+\tau}{\tau}*\pi) \f}
 * 
 * @code
 *             case 2:{
 * 
 *                 noise = (gen() * S);
 *                 t0 = P2*double(col);
 * 
 *                 for(int k = 0; k<col; k++){
 * 
 *                     t = 12;
 *                     if(k < t0-t/2 || k > t0+t/2)
 *                         waveforms(j,k) = 0.0;
 *                     else
 *                         waveforms(j,k) = (1 + noise) * std::pow(cos(((k - t0 + t * noise) / t) * M_PI), 2);
 *                     waveforms(j,k)   += gen() * S;
 *                 }
 * 
 * 
 *                 break;
 *             }
 * @endcode
 * 
 * \f{equation} Y = \frac{1}{1+\left(\frac{x-x_{0}}{\sigma_{0}}\right)^2} - \frac{1}{1+\left(\frac{x-x_{1}}{\sigma_{1}}\right)^2}    + R \f}
 * 
 * @code
 *             case 3:{
 * 
 *                 noise = 1+(gen()*S)*0.5;
 * 
 *                 X1 = P3*noise;
 *                 sig1 = P5*noise;
 *                 t = noise;
 *                 X2 = P4*noise;
 *                 sig2 = P6*noise;
 *                 l = 0.;
 * 
 *                 for(int k = 0; k<col; k++){
 *                     Y = gen()*S;
 *                     l = k/double(col);
 *                     waveforms(j,k) = (1+Y)+(1/t)*(((1/(1+std::pow(((l-X1)/sig1),2))) - (1/(1+std::pow(((l-X2)/sig2),2)))));
 *                 }
 *                 break;
 *             }
 *             default:{
 *                 std::cout<<"illegal waveform choosen"<<std::endl;
 *             }
 * 
 *             }
 *         }
 *     }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionnormalize_waveforms"></a> 
 * <h4>Function: normalize_waveforms</h4>
 * This function normalizes the generated waveforms in respect to the X axis
 * 
 * 
 * @param n_rows : number of rows/datapoints.
 * @param n_cols : number of columns (space dimension).
 * @param waveforms : the matrix containing the generated waveforms.
 * 
 * @code
 *     template <typename D>
 *     static void normalize_waveforms(unsigned int n_rows, unsigned int n_cols,
 *                                     dealii::FullMatrix<D> &waveforms)
 *     {
 *         double row_sum, row_mean;
 *         for(int i = 0; i < n_rows; i++){
 *             row_sum = 0., row_mean = 1.;
 *             for(int j = 0; j < n_cols; j++){
 *                 row_sum += waveforms(i,j);
 *             }
 *             row_mean = row_sum/(double)n_cols;
 *             for(int j = 0; j < n_cols; j++){
 *                 waveforms(i,j) -= row_mean;
 *             }
 *         }
 *     }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionoriginal_waveforms"></a> 
 * <h4>Function: original_waveforms</h4>
 * This function generates one waveform from each one of the three available forms
 * those will be used as the real means for the ShohamEM algorithm
 * 
 * 
 * @param col : number of columns (space dimension).
 * @param n_forms : number of forms choosen by the user.
 * @param check_waves : the matrix where the data is written to.
 * @param waves_parameter : conataints the parameter needed for generating the diffrent waveforms.
 * @param n_waves : number of waves from each waveform.
 * 
 * @code
 *     template <typename D>
 *     static void original_waveforms(unsigned int col, int n_forms,
 *                                    dealii::FullMatrix<D> &check_waves,
 *                                    const std::vector<double> &waves_parameter,
 *                                    const std::vector<int> &n_waves)
 *     {
 *         double X, X1, X2, sig, sig1, sig2, nu, t, noise, l;
 *         double P0,P1,P2,P3,P4,P5,P6;
 * 
 *         P0 = waves_parameter[0];
 *         P1 = waves_parameter[1];
 *         P2 = waves_parameter[2];
 *         P3 = waves_parameter[3];
 *         P4 = waves_parameter[4];
 *         P5 = waves_parameter[5];
 *         P6 = waves_parameter[6];
 * 
 *         int w1, w2, w3;
 * 
 *         w1 = n_waves[0];
 *         w2 = n_waves[1];
 *         w3 = n_waves[2];
 * 
 *         for(int i=0; i<n_forms;i++){
 *             if(w1!= 0 && i<n_forms){
 *                 X = P0;
 *                 sig = P1;
 *                 t = 0.1;
 *                 l = 0.;
 * 
 *                 for(int k = 0; k<col; k++){
 *                     l = k/double(col);
 *                     check_waves(i,k) = (1/t)*((l-X)*exp(-(std::pow(((l-X)/(sig*sqrt(2))),2))));
 *                 }
 *                 w1 = 0;
 *                 i ++;
 *             }
 *             if(w2!=0 && i<n_forms){
 * 
 *                 nu = P2;
 *                 X = (1/(nu*2));
 *                 l= 0.;
 * 
 *                 for(int k = 0; k<col; k++){
 *                     l = k/double(col) + noise;
 *                     if(l > 2*X|| l < 0 ){
 *                         check_waves(i,k) = 0.;
 *                     }else{
 *                         check_waves(i,k) = std::pow(sin(nu*M_PI*l),2);
 *                     }
 *                 }
 *                 w2 =0;
 *                 i++;
 *             }
 *             if(w3!=0 && i<n_forms){
 * 
 *                 X1 = P3;
 *                 sig1 = P5;
 *                 X2 = P4;
 *                 sig2 = P6;
 *                 l = 0.;
 * 
 *                 for(int k = 0; k<col; k++){
 *                     l = k/double(col);
 *                     check_waves(i,k) = (((1/(1+std::pow(((l-X1)/sig1),2))) - (1/(1+std::pow(((l-X2)/sig2),2)))));
 *                 }
 *                 w3 = 0;
 *                 i++;
 *             }
 *         }
 *     }
 * };
 * 
 * 
 * #endif // WAVEFORMS_GENERATOR_H
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef CUDA_KERNEL_STEP_7_CU_H
 * #define CUDA_KERNEL_STEP_7_CU_H
 * 
 * @endcode
 * 
 * 
 * <a name="Declarationsofthekernelwrapperfunctions"></a> 
 * <h3>Declarations of the kernel-wrapper-functions</h3>
 * 
 * @code
 * template <typename T>
 * struct Kernel {
 * 
 *     void subtract_array_from_matrix(int n, int p, int j, int g,
 *                                     const T * X,
 *                                     const T * mu_j,
 *                                     T * M);
 * 
 *     void MM_scalar(int n, int p,
 *                    const T * A,
 *                    const T * B,
 *                    T * result);
 * 
 *     T Det(int p, int MS, const T * A);
 * 
 * 
 * 
 *     void assemble_P_ij(T gamma,
 *                        T v,
 *                        int p,int n, int g,
 *                        const T * Delta,
 *                        const T * delta,
 *                        T * P_ij);
 * 
 * 
 *     void assemble_u_ij(T v,
 *                        int p,int n, int g,
 *                        const T * delta,
 *                        T * u_ij);
 * 
 *     void assemble_z_ij(int n, int g,
 *                        const T * P_ij,
 *                        const T * c,
 *                        const T * pi,
 *                        T * z_ij);
 * 
 *     void reduce(unsigned int n,
 *                 const T * c,
 *                 T * result,
 *                 int threads,
 *                 int blocks);
 * 
 *     double calculate_y(int n,
 *                        int g,
 *                        T v,
 *                        T digamma,
 *                        const T * z_ij,
 *                        const T * delta_ij,
 *                        const T * u_ij);
 * 
 *     double sum_array(int n,
 *                      const T * input,
 *                      bool log);
 * 
 *     void sum_matrix_cols(int n, int g,
 *                          const T * input,
 *                          T * output,
 *                          bool log);
 * 
 *     double sum_array(int n,
 *                      int g,
 *                      const T * input,
 *                      bool log);
 * 
 *     void sum_array2D(int n,
 *                      int g,
 *                      const T * z_ij,
 *                      const T * u_ij,
 *                      T * output);
 * 
 *     void MxM(int n,int g,
 *              const T * z_ij,
 *              const T * u_ij,
 *              T * q_ij);
 * 
 *     void SMX(int n, int p,
 *              const T * q_ij,
 *              const T C_j,
 *              const int j,
 *              const T * X,
 *              const T * mu_j,
 *              T * M);
 * 
 * };
 * 
 * #endif // CUDA_KERNEL_STEP_7_CU_H
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * @endcode
 * 
 * Header-File of the CUDA utility-Library.
 * 
 * 
 * 
 * @code
 * #include <step-7/cuda_kernel_wrapper_step-7.cu.h>
 * 
 * #include <math.h>
 * #include <stdio.h>
 * #include <iostream>
 * #include <fstream>
 * 
 * float Giga_bytes = std::pow(1024.,3);
 * 
 * cudaEvent_t beginEvent;
 * cudaEvent_t endEvent;
 * 
 * std::ofstream dt("output/Memory_Bandwidth.txt");
 * 
 * @endcode
 * 
 * 
 * <a name="Kernel__subtract_array_from_matrix"></a> 
 * <h3>Kernel: __subtract_array_from_matrix</h3>
 * Subtracts the array @p mu_j from each row of the matrix @p X .
 * The results are saved in the output matrix @p M .
 * 
 * @code
 * template <typename T>
 * __global__ void __subtract_array_from_matrix(int n, int p, int j, int g,
 *                                              const T * X,
 *                                              const T * mu_j,
 *                                              T * M)
 * {
 *     int r = blockIdx.y * blockDim.x + threadIdx.x;
 *     int c = blockIdx.x;
 * 
 *     if(r < n && c < p){
 *         M[c*n+r] = X[c*n+r]-mu_j[j+c*g];
 *     }
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Kernel__MM_scalar"></a> 
 * <h3>Kernel: __MM_scalar</h3>
 * calculates the Matrix Matrix Scalarproduct by calculating the Scalarproduct of each row from the first matrix
 * with it corresponding column from the second matrix
 * 
 * @code
 * template <typename T>
 * __global__ void __MM_scalar(int n, int p,
 *                             const T * A,
 *                             const T * B,
 *                             T * result)
 * {
 *     int r = blockIdx.x * blockDim.x + threadIdx.x;
 * 
 *     T sum = 0.;
 * 
 *     if(r < n){
 *         for(int j= 0; j<p; j++){
 *             sum += A[j*n+r]*B[r*p+j];
 *         }
 *         result[r] = sum;
 *     }
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Kernel__Det"></a> 
 * <h3>Kernel: __Det</h3>
 * calculate the determinate
 * 
 * @code
 * template <typename T>
 * __global__ void __Det(int p, int MS,
 *                       const T * A,
 *                       T * result)
 * {
 *     T det = 1.;
 * 
 *     for(int i = 0; i<p; i++){
 *         det *= A[i*MS+i];
 *     }
 * 
 *     *result = det;
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Kernel__P_ij"></a> 
 * <h3>Kernel: __P_ij</h3>
 * calculates the @p P_ij matrix by solving the follwing equation:
 * \f{equation} P_{ij} = \frac{\gamma}{\Delta_j}\frac{1}{(1+\frac{\delta_{ij}}{\nu})^{(\nu+p)/2}}\f}
 * and \f{equation} \gamma = \frac{\Gamma(\frac{\nu+p}{2})}{\Gamma(\frac{\nu}{2})(\pi\nu)^{p/2}}\f}
 * 
 * @code
 * template <typename T>
 * __global__ void __P_ij(T gamma,
 *                        T v,
 *                        int p, int n, int g,
 *                        const T * Delta,
 *                        const T * delta,
 *                        T * P_ij)
 * {
 *     int i = blockIdx.y * blockDim.x + threadIdx.x;
 *     int j = blockIdx.x;
 * 
 * 
 *     if(i < n && j < g){
 *         P_ij[j*n+i] = gamma/Delta[j]/(pow((1.+delta[j*n+i]/v),(v+p)/2.));
 *     }
 * 
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Kernel__u_ij"></a> 
 * <h3>Kernel: __u_ij</h3>
 * calculates the @p u_ij matrix by solving the following equation:
 * \f{equation} u_{ij} = \frac{p+2}{\delta_{ij}+\nu}\f}
 * 
 * @code
 * template <typename T>
 * __global__ void __u_ij(T v,
 *                        int p, int n, int g,
 *                        const T * delta,
 *                        T * u_ij)
 * {
 *     int i = blockIdx.y * blockDim.x + threadIdx.x;
 *     int j = blockIdx.x;
 * 
 * 
 *     if(i < n && j < g){
 *         u_ij[j*n+i] = (p+2)/(delta[j*n+i]+v);
 *     }
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Kernel__z_ij"></a> 
 * <h3>Kernel: __z_ij</h3>
 * calculates the @p z_ij matrix by solving the following equation:
 * \f{equation} u_{ij} = \frac{P_{ij}\pi_j}{c}\f}
 * and \f{equation} c = \sum_{l=1}^{g}P_{il}\pi_l\f}
 * 
 * @code
 * template <typename T>
 * __global__ void __z_ij(int n, int g,
 *                        const T * P_ij,
 *                        const T * c,
 *                        const T * pi,
 *                        T * z_ij)
 * {
 *     __shared__ T cs[512];
 *     __shared__ T pi_s[32];
 * 
 *     int i = threadIdx.x + blockDim.x * blockIdx.x;
 *     int j = threadIdx.y;
 * 
 * 
 *     if (j == 0){
 *         cs[threadIdx.x] = c[i];
 *     }
 * 
 * @endcode
 * 
 * using the threadIdx.x to load data into the sharedmemory because it runs faster than threadIdx.y
 * 
 * @code
 *     if(threadIdx.x < g && j == 0){
 *         pi_s[threadIdx.x] = pi[threadIdx.x];
 *     }
 * 
 *     __syncthreads();
 * 
 *     if(i < n && j < g){
 *         z_ij[j*n+i] = (P_ij[j*n+i]*pi_s[j])/cs[threadIdx.x];
 *     }
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionreduction"></a> 
 * <h4>Function: reduction</h4>
 * 
 * 
 * performs the last reduction on the sharedmemory <br>
 * if thread_block isn't fill_fledged, it could cause some problems
 * 
 * @code
 * template<typename T>
 * __device__ void reduction(volatile T * sum, int tid)
 * {
 *     for(int k = blockDim.x/2; k>0; k/=2){
 *         if(tid < k){
 *             sum[tid] += sum[tid + k];
 *         }
 *         __syncthreads();
 *     }
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Kernel__sum_array"></a> 
 * <h3>Kernel: __sum_array</h3>
 * 
 * 
 * sums up the entries of the input array in @p result
 * 
 * @code
 * template <typename T>
 * __global__ void __sum_array(int n,
 *                             const T * input,
 *                             T * result)
 * {
 *     __shared__ T sum[512];
 * 
 *     sum[threadIdx.x] = 0.;
 *     __syncthreads();
 * 
 * 
 *     int i = threadIdx.x + blockDim.x * blockIdx.x;
 * 
 *     if(i>=n)
 *         return;
 * 
 *     sum[threadIdx.x] += input[i];
 *     __syncthreads();
 * 
 *     reduction(sum, threadIdx.x);
 * 
 *     if(threadIdx.x == 0){
 *         result[blockIdx.x]= sum[0];
 *     }
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Kernel__sum_log_array"></a> 
 * <h3>Kernel: __sum_log_array</h3>
 * 
 * 
 * sums up the logarithm of the entries of the input array in @p result
 * 
 * @code
 * template <typename T>
 * __global__ void __sum_log_array(int n,
 *                                 const T * input,
 *                                 T * result)
 * {
 *     __shared__ T sum[512];
 * 
 *     sum[threadIdx.x] = 0.;
 *     __syncthreads();
 * 
 * 
 *     int i = threadIdx.x + blockDim.x * blockIdx.x;
 * 
 *     if(i>=n)
 *         return;
 * 
 *     sum[threadIdx.x] += log(input[i]);
 *     __syncthreads();
 * 
 *     reduction(sum, threadIdx.x);
 * 
 *     if(threadIdx.x == 0){
 *         result[blockIdx.x]= sum[0];
 *     }
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Kernel__sum_array2D"></a> 
 * <h3>Kernel: __sum_array2D</h3>
 * 
 * 
 * calculates the row sum of the multiplication of @p z_ij with @p u_ij
 * 
 * @code
 * template <typename T>
 * __global__ void __sum_array2D(int n, int g,
 *                               const T * z_ij,
 *                               const T * u_ij,
 *                               T * output)
 * {
 *     __shared__ T sum[512];
 * 
 *     sum[threadIdx.x] = 0.;
 *     __syncthreads();
 * 
 * 
 *     int i = threadIdx.x + blockDim.x * blockIdx.x;
 *     int j = blockIdx.y;
 * 
 *     if(i>=n)
 *         return;
 * 
 *     int global_index = j*n+i;
 *     sum[threadIdx.x] += z_ij[global_index]* u_ij[global_index];
 *     __syncthreads();
 * 
 *     reduction(sum, threadIdx.x);
 * 
 *     if(threadIdx.x == 0){
 *         output[blockIdx.y*gridDim.x+blockIdx.x]= sum[0];
 *     }
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Kernel__y"></a> 
 * <h3>Kernel: __y</h3>
 * 
 * 
 * calculates y with:
 * \f{equation} y = \sum_{i=1}^{n}\sum_{j=1}^{g}z_{ij}[\psi+log(\frac{2}{\delta_{ij}+\nu_{old}})-u_{ij}]\f}
 * and \f{equation} \psi = \psi(\frac{p+\nu{old}}{2}) \f} with $\psi$ = the digamma function.
 * 
 * @code
 * template <typename T>
 * __global__ void __y(int n, int g, T v, T digamma,
 *                     const T * z_ij,
 *                     const T * delta_ij,
 *                     const T * u_ij,
 *                     T * results)
 * {
 *     __shared__ T sum[256];
 * 
 *     sum[threadIdx.x] = 0.;
 *     __syncthreads();
 * 
 * 
 *     int i = threadIdx.x + blockDim.x * blockIdx.x;
 * 
 *     if(i>=n)
 *         return;
 * 
 *     for(int j = 0; j<g; j++){
 *         int global_index = j*n+i;
 *         sum[threadIdx.x] += z_ij[global_index]*(digamma+log(2/(delta_ij[global_index]+v))-u_ij[global_index]);
 *     }
 *     __syncthreads();
 * 
 *     reduction(sum, threadIdx.x);
 * 
 *     if(threadIdx.x == 0){
 *         results[blockIdx.x]= sum[0];
 *     }
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Kernel__MxM"></a> 
 * <h3>Kernel: __MxM</h3>
 * calculates the prduct of each entry of the first matrix with the corresponding one inn the second matrix.
 * the results are saved in the new matrix @p q_ij.
 * 
 * @code
 * template <typename T>
 * __global__ void __MxM(int n, int g,
 *                       const T * z_ij,
 *                       const T * u_ij,
 *                       T * q_ij)
 * {
 *     int i = blockIdx.y * blockDim.x + threadIdx.x;
 *     int j = blockIdx.x;
 * 
 * 
 *     if(i < n && j < g){
 *         q_ij[j*n+i] = z_ij[j*n+i]*u_ij[j*n+i];
 *     }
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Kernel__SMX"></a> 
 * <h3>Kernel: __SMX</h3>
 * calculates the following:
 * \f{equation} M = \sum_{i=1}^{n}\sqrt{\frac{q_{ij}}{C_j}}(\vec{x_i}-\vec{\mu_j})\f}
 * with \f{equation} q_{ij} = z_{ij}u_{ij} \f}
 * and @p C_j = the j_th array from the matrix C which is
 * \f{equation} C = \sum_{i=1}^{n}z_{ij}u_{ij}\f}
 * therefore by calculating $M^T * M$ we get
 * \f{equation} \Sigma_j = \sum_{i=1}^{n}\frac{\left(z_{ij}u_{ij}\right)(\vec{x_i}-\vec{\mu_j})(\vec{x_i}-\vec{\mu_j})^T}{\sum_{i=1}^{n}z_{ij}u_{ij}}\f}
 * 
 * @code
 * template <typename T>
 * __global__ void __SMX(int n, int p,
 *                       const T * q_ij,
 *                       const T inv_C_j,
 *                       const int j,
 *                       const T * X,
 *                       const T * mu_j,
 *                       T * M)
 * {
 *     int r = blockIdx.y * blockDim.x + threadIdx.x;
 *     int c = blockIdx.x;
 * 
 * 
 *     if(r < n && c < p){
 *         M[c*n+r] = (X[c*n+r]-mu_j[c]) * sqrt(q_ij[j*n+r] * inv_C_j);
 *     }
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionsubtract_array_from_matrix"></a> 
 * <h4>Function:subtract_array_from_matrix</h4>
 * 
 * 
 * Wrapper-Function for the Kernel. Sets the threads distribution and runs them.
 * @param n : number of rows in the input Matrix @p X (represents the number of datapoints).
 * @param p : number of rows (represents the space dimension).
 * @param j : index of the array @p mu_j.
 * @param g : number of rows in the array  @p mu_j (represents the number of cluster centers).
 * @param X : input matrix with @p n rows and @p p columns.
 * @param mu_j : the array which is subtracted from the matrix @p X .
 * @param M : output result matrix.
 * 
 * @code
 * template<typename T>
 * void Kernel<T>::subtract_array_from_matrix(int n, int p, int j, int g,
 *                                            const T * X,
 *                                            const T * mu_j,
 *                                            T * M)
 * {
 *     int block_size = 512;
 *     int gridDim_y = (n+block_size-1)/block_size;
 * 
 *     dim3 num_blocks(p, gridDim_y);
 * 
 *     double time_sum = 0;
 *     int n_operations = 1;
 *     int n_memory_access = 3;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * 
 *     dt <<"Kernel name"
 *       <<"\t"
 *      <<"Time in Milliseconds"
 *     <<"\t\t"
 *     <<"bandwidth in GB/s"
 *     <<"\t\t"
 *     <<"Gflops/Second"
 *     << std::endl;
 *     dt << "------------------------------------------------------------------------------------------------------------"
 *        <<std::endl;
 * 
 * @endcode
 * 
 * printf("Kernel Running ");
 * 
 * @code
 *     for(int s = 0; s < N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         __subtract_array_from_matrix<<<num_blocks, block_size>>>(n, p, j, g,
 *                                                                  X,
 *                                                                  mu_j,
 *                                                                  M);
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * @endcode
 * 
 * if((s%10) == 0) printf(". ");
 * 
 * @code
 *     }
 * @endcode
 * 
 * printf("\n");
 * 
 * @code
 *     time_sum /= N_REPS;
 *     dt << "subtract_kernel:::::"
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*p*sizeof(T))*n_memory_access)/(time_sum/1000))/ Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*g*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="FunctionMM_scalar"></a> 
 * <h4>Function: MM_scalar</h4>
 * 
 * 
 * Wrapper-Function for the Kernel. Sets the threads distribution and runs them.
 * 
 * 
 * @param n : number of rows in the input matrices.
 * @param p : number of columns in the input matrices (represents the space dimension).
 * @param A : first input matrix of the size @p n X @p p.
 * @param B : second input matrix.
 * @param output : result vector.
 * 
 * @code
 * template<typename T>
 * void Kernel<T>::MM_scalar(int n, int p,
 *                           const T * A,
 *                           const T * B,
 *                           T * output)
 * {
 *     int block_size = 512;
 *     int num_blocks = (n+block_size-1)/block_size;
 * 
 *     double time_sum = 0;
 *     long n_operations = 2*p;
 *     int n_memory_access = 2*p+1;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * @endcode
 * 
 * printf("Kernel Running ");
 * 
 * @code
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 *         __MM_scalar<<<num_blocks, block_size>>>(n, p,
 *                                                 A,
 *                                                 B,
 *                                                 output);
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * @endcode
 * 
 * if((s%10) == 0) printf(". ");
 * 
 * @code
 *     }
 * @endcode
 * 
 * printf("\n");
 * 

 * 
 * 
 * @code
 *     time_sum /= N_REPS;
 *     dt << "MM_scalar:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*p*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*p)*n_operations)/(time_sum/1000)/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="FunctionDet"></a> 
 * <h4>Function: Det</h4>
 * 
 * 
 * Wrapper-Function for the Kernel. Sets the threads distribution and runs them.
 * 
 * 
 * @param p : size of the matrixview @p A .
 * @param MS : The size of the original matrix where @p A is stored.
 * @param A : input matrixview on the original matrix.
 * @return det : returns the determinant of the matrix @p A .
 * 
 * @code
 * template<typename T>
 * T
 * Kernel<T>::Det(int p, int MS, const T * A)
 * {
 *     int block_size = 1;
 *     int num_blocks = 1;
 * 
 *     double time_sum = 0;
 *     int n_operations = p;
 *     int n_memory_access = p+1;
 * 
 * 
 *     T det;
 *     T * det_d;
 *     cudaMalloc(&det_d, sizeof(T));
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * @endcode
 * 
 * printf("Kernel Running ");
 * 
 * @code
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         __Det<<<num_blocks, block_size>>>(p, MS, A, det_d);
 * 
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * @endcode
 * 
 * if((s%10) == 0) printf(". ");
 * 

 * 
 * 
 * @code
 *     }
 * @endcode
 * 
 * printf("\n");
 * 
 * @code
 *     time_sum /= N_REPS;
 *     dt << "Det_Kernel:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((p*p*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((p*p*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 *     cudaThreadSynchronize();
 * 
 *     cudaMemcpy(&det, det_d,
 *                sizeof(T), cudaMemcpyDeviceToHost);
 * 
 *     cudaFree(det_d);
 * 
 *     return det;
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionassemble_P_ij"></a> 
 * <h4>Function: assemble_P_ij</h4>
 * 
 * 
 * Wrapper-Function for the Kernel. Sets the threads distribution and runs them.
 * 
 * 
 * @param gamma : value of gamma calculated on the CPU.
 * @param v : DOF parameter.
 * @param p : number of columns.
 * @param n : number of rows.
 * @param g : number of cluster centers.
 * @param Delta : input $Delta_j$ array.
 * @param delta : input $delta_{ij}$ matrix containing the Mahalanobis squared distance betweet each datapoint the the real cluster centers.
 * @param P_ij : result matrix.
 * 
 * @code
 * template<typename T>
 * void
 * Kernel<T>::assemble_P_ij(T gamma,
 *                          T v,
 *                          int p, int n, int g,
 *                          const T * Delta,
 *                          const T * delta,
 *                          T * P_ij)
 * {
 *     int block_size = 512;
 *     int gridDim_y = (n+block_size-1)/block_size;
 * 
 *     dim3 num_blocks(g, gridDim_y);
 * 
 * 
 *     double time_sum = 0;
 *     int n_operations = 4;
 *     int n_memory_access = 3;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * @endcode
 * 
 * printf("Kernel Running ");
 * 
 * @code
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 * 
 *         __P_ij<<<num_blocks, block_size>>>(gamma,
 *                                            v,
 *                                            p, n, g,
 *                                            Delta,
 *                                            delta,
 *                                            P_ij);
 * 
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * @endcode
 * 
 * if((s%10) == 0) printf(". ");
 * 

 * 
 * 
 * @code
 *     }
 * @endcode
 * 
 * printf("\n");
 * 
 * @code
 *     time_sum /= N_REPS;
 *     dt << "assemble_P_ij:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*g*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*g*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionassemble_u_ij"></a> 
 * <h4>Function: assemble_u_ij</h4>
 * 
 * 
 * Wrapper-Function for the Kernel. Sets the threads distribution and runs them.
 * 
 * 
 * @param v : DOF parameter.
 * @param p : number of columns.
 * @param n : number of rows.
 * @param g : number of cluster centers.
 * @param delta : input Mahalanobis squared distances matrix $\delta_{ij}$.
 * @param u_ij : output weights matrix.
 * 
 * @code
 * template<typename T>
 * void
 * Kernel<T>::assemble_u_ij(T v,
 *                          int p, int n, int g,
 *                          const T * delta,
 *                          T * u_ij)
 * {
 *     int block_size = 512;
 *     int gridDim_y = (n+block_size-1)/block_size;
 * 
 *     dim3 num_blocks(g, gridDim_y);
 * 
 * 
 *     double time_sum = 0;
 *     int n_operations = 2;
 *     int n_memory_access = 2;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * @endcode
 * 
 * printf("Kernel Running ");
 * 

 * 
 * 
 * @code
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         __u_ij<<<num_blocks, block_size>>>(v,
 *                                            p, n, g,
 *                                            delta,
 *                                            u_ij);
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * @endcode
 * 
 * if((s%10) == 0) printf(". ");
 * 

 * 
 * 
 * @code
 *     }
 * @endcode
 * 
 * printf("\n");
 * 
 * @code
 *     time_sum /= N_REPS;
 *     dt << "assemble_u_ij:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*g*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*g*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionassemble_z_ij"></a> 
 * <h4>Function: assemble_z_ij</h4>
 * 
 * 
 * Wrapper-Function for the Kernel. Sets the threads distribution and runs them.
 * 
 * 
 * @param n : number of rows in the @p P_ij matrix and c array
 * @param g : number of cluster centers (in this case the number of rows in the $\pi$ array
 * @param P_ij : input matrix representing the distribution of the jth spike.
 * @param c : c array.
 * @param pi : proportion array $\pi$
 * @param z_ij : output memberships matrix.
 * 
 * @code
 * template<typename T>
 * void
 * Kernel<T>::assemble_z_ij(int n, int g,
 *                          const T * P_ij,
 *                          const T * c,
 *                          const T * pi,
 *                          T * z_ij)
 * {
 *     int blockDim_x = max(16/g, 1) * 32;
 *     dim3 block_size(blockDim_x,g);
 *     int num_blocks = (n+blockDim_x-1)/blockDim_x;
 * 
 * 
 *     double time_sum = 0;
 *     int n_operations = 2;
 *     int n_memory_access = 2;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * @endcode
 * 
 * printf("Kernel Running ");
 * 

 * 
 * 
 * @code
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 *         __z_ij<<<num_blocks, block_size>>>(n, g,
 *                                            P_ij,
 *                                            c,
 *                                            pi,
 *                                            z_ij);
 * 
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * @endcode
 * 
 * if((s%10) == 0) printf(". ");
 * 

 * 
 * 
 * @code
 *     }
 * @endcode
 * 
 * printf("\n");
 * 
 * @code
 *     time_sum /= N_REPS;
 *     dt << "assemble_z_ij:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << ((((n*g*sizeof(T))*n_memory_access)+n+g)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*g*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functioncalculate_y"></a> 
 * <h4>Function: calculate_y</h4>
 * 
 * 
 * Wrapper-Function for the Kernel. Sets the threads distribution and runs them.
 * @param n : number of rows (number of datapoints).
 * @param g : number of columns (in this case the number of cluster centers).
 * @param v : DOF parameter.
 * @param digamma : the digamma function calculated in the CPU.
 * @param z_ij : input memberships matrix.
 * @param delta_ij : input Mahalanobis squared distances matrix.
 * @param u_ij : input weight matrix.
 * @return y : the value y.
 * 
 * @code
 * template<typename T>
 * double
 * Kernel<T>::calculate_y(int n, int g, T v, T digamma,
 *                        const T * z_ij,
 *                        const T * delta_ij,
 *                        const T * u_ij)
 * {
 *     int block_size = 256;
 *     int num_blocks = (n+block_size-1)/block_size;
 * 
 *     T y;
 * 
 *     T * results_h = new T[num_blocks];
 *     T * results_d;
 *     cudaMalloc(&results_d, num_blocks*sizeof(T));
 * 
 *     double time_sum = 0;
 *     int n_operations = 6*g + sqrt(block_size/2)+1;
 *     int n_memory_access = 4;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * @endcode
 * 
 * printf("Kernel Running ");
 * 

 * 
 * 
 * @code
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         __y<<<num_blocks, block_size>>>(n, g, v, digamma,
 *                                         z_ij,
 *                                         delta_ij,
 *                                         u_ij,
 *                                         results_d);
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * @endcode
 * 
 * if((s%10) == 0) printf(". ");
 * 

 * 
 * 
 * @code
 *     }
 * @endcode
 * 
 * printf("\n");
 * 
 * @code
 *     time_sum /= N_REPS;
 *     dt << "calculate_y:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*g*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*g*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 *     cudaThreadSynchronize();
 * 
 *     cudaMemcpy(results_h, results_d,
 *                num_blocks*sizeof(T), cudaMemcpyDeviceToHost);
 * 
 *     for(int i=0; i<num_blocks; i++){
 *         y += results_h[i];
 *     }
 * 
 *     return y;
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionsum_array"></a> 
 * <h4>Function: sum_array</h4>
 * 
 * 
 * Wrapper-Function for the Kernel. Sets the threads distribution and runs them.
 * @param n : number of elements in the input array.
 * @param input : input array with @p n entries.
 * @param log : determin summing up the entries or the logarithm of the entries.
 * @return result : the sum of the input entries.
 * 
 * @code
 * template<typename T>
 * double
 * Kernel<T>::sum_array(int n,
 *                      const T * input,
 *                      bool log)
 * {
 *     double result = 0.;
 *     int block_size = 512;
 *     int num_blocks = (n+block_size-1)/block_size;
 * 
 * 
 *     T * results_h = new T[num_blocks];
 *     T * results_d;
 *     cudaMalloc(&results_d, num_blocks*sizeof(T));
 * 
 *     double time_sum = 0;
 *     int n_operations = sqrt(block_size/2)+1+1;
 *     int n_memory_access = 2;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * @endcode
 * 
 * printf("Kernel Running ");
 * 

 * 
 * 
 * @code
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         if(log){
 *             __sum_log_array<<<num_blocks, block_size>>>(n,input,results_d);
 *         }else{
 *             __sum_array<<<num_blocks, block_size>>>(n,input,results_d);
 *         }
 * 
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * @endcode
 * 
 * if((s%10) == 0) printf(". ");
 * 

 * 
 * 
 * @code
 *     }
 * @endcode
 * 
 * printf("\n");
 * 
 * @code
 *     time_sum /= N_REPS;
 *     dt << "sum_array:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 *     cudaThreadSynchronize();
 * 
 *     cudaMemcpy(results_h, results_d,
 *                num_blocks*sizeof(T), cudaMemcpyDeviceToHost);
 * 
 * @endcode
 * 
 * sums up the results of each block to one value
 * 
 * @code
 *     for(int i=0; i<num_blocks; i++){
 *         result += results_h[i];
 *     }
 * 
 *     return result;
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionsum_matrix_cols"></a> 
 * <h4>Function: sum_matrix_cols</h4>
 * 
 * 
 * Wrapper-Function for the Kernel. Sets the threads distribution and runs them.
 * 
 * 
 * @param n : number of rows.
 * @param g : number of cluster centers (columns).
 * @param input : input matrix with n rows and g columns.
 * @param output : output array containing the results of the sum of each matrix column.
 * @param log : determin summing up the logarithm of the matrix entries or not.
 * 
 * @code
 * template<typename T>
 * void
 * Kernel<T>::sum_matrix_cols(int n, int g,
 *                            const T * input,
 *                            T * output,
 *                            bool log)
 * {
 *     int block_size = 512;
 *     int num_blocks = (n+block_size-1)/block_size;
 * 
 * 
 *     T * results_h = new T[num_blocks*g];
 *     T * results_d;
 *     cudaMalloc(&results_d, num_blocks*g*sizeof(T));
 * 
 * @endcode
 * 
 * the Kernel is called g times for each row of the matrix
 * 
 * @code
 *     for(int j = 0;j < g; j++){
 *         if(log){
 *             __sum_log_array<<<num_blocks, block_size>>>(n,input+(n*j),results_d+(j*num_blocks));
 *         }else{
 *             __sum_array<<<num_blocks, block_size>>>(n,input+(n*j),results_d+(j*num_blocks));
 *         }
 *     }
 *     cudaThreadSynchronize();
 * 
 *     cudaMemcpy(results_h, results_d,
 *                num_blocks*g*sizeof(T), cudaMemcpyDeviceToHost);
 * 
 * @endcode
 * 
 * for each cluster center sums up the results of each block to one value.
 * 
 * @code
 *     for(int j = 0; j<g; j++){
 *         for(int i=0; i<num_blocks; i++){
 *             output[j] += results_h[i+j*num_blocks];
 *         }
 *     }
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionsum_array2D"></a> 
 * <h4>Function: sum_array2D</h4>
 * 
 * 
 * Wrapper-Function for the Kernel. Sets the threads distribution and runs them.
 * 
 * 
 * @param n : number of rows.
 * @param g : number of columns.
 * @param z_ij : input memberships matrix.
 * @param u_ij : input weights matrix.
 * @param output : output array withthe sum of each column in each entry.
 * 
 * @code
 * template<typename T>
 * void
 * Kernel<T>::sum_array2D(int n, int g,
 *                        const T * z_ij,
 *                        const T * u_ij,
 *                        T * output)
 * {
 *     int block_size = 512;
 *     int gridDim_x = ((n+block_size-1)/block_size);
 *     dim3 num_blocks (gridDim_x, g);
 * 
 * 
 *     T * results_h = new T[gridDim_x*g];
 *     T * results_d;
 *     cudaMalloc(&results_d, gridDim_x*g*sizeof(T));
 * 
 *     double time_sum = 0;
 *     int n_operations = 2 + sqrt(block_size/2)+1;
 *     int n_memory_access = 2;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * @endcode
 * 
 * printf("Kernel Running ");
 * 

 * 
 * 
 * @code
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         __sum_array2D<<<num_blocks, block_size>>>(n,g,z_ij, u_ij,results_d);
 * 
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * @endcode
 * 
 * if((s%10) == 0) printf(". ");
 * 

 * 
 * 
 * @code
 *     }
 * @endcode
 * 
 * printf("\n");
 * 
 * @code
 *     time_sum /= N_REPS;
 *     dt << "sum_array2D:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*g*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*g*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 *     cudaThreadSynchronize();
 * 
 *     cudaMemcpy(results_h, results_d,
 *                gridDim_x*g*sizeof(T), cudaMemcpyDeviceToHost);
 * 
 *     for(int i=0; i<g; i++){
 *         for(int j= 0; j<gridDim_x; j++){
 *             output[i] += results_h[i*gridDim_x+j];
 *         }
 *     }
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="FunctionMxM"></a> 
 * <h4>Function: MxM</h4>
 * 
 * 
 * Wrapper-Function for the Kernel. Sets the threads distribution and runs them.
 * 
 * 
 * @param n : number of rows.
 * @param g : number of columns.
 * @param z_ij : first input matrix.
 * @param u_ij : second input matrix.
 * @param q_ij : output results matrix.
 * 
 * @code
 * template<typename T>
 * void
 * Kernel<T>::MxM(int n, int g,
 *                const T * z_ij,
 *                const T * u_ij,
 *                T * q_ij)
 * {
 *     int block_size = 512;
 *     int gridDim_y = (n+block_size-1)/block_size;
 * 
 *     dim3 num_blocks(g, gridDim_y);
 * 
 * 
 *     double time_sum = 0;
 *     int n_operations = 1;
 *     int n_memory_access = 3;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * @endcode
 * 
 * printf("Kernel Running ");
 * 

 * 
 * 
 * @code
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         __MxM<<<num_blocks, block_size>>>(n, g,
 *                                           z_ij,
 *                                           u_ij,
 *                                           q_ij);
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * @endcode
 * 
 * if((s%10) == 0) printf(". ");
 * 

 * 
 * 
 * @code
 *     }
 * @endcode
 * 
 * printf("\n");
 * 
 * @code
 *     time_sum /= N_REPS;
 *     dt << "MxM_Kernel:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*g*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*g*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="FunctionSMX"></a> 
 * <h4>Function: SMX</h4>
 * 
 * 
 * Wrapper-Function for the Kernel. Sets the threads distribution and runs them.
 * @param n : number of rows.
 * @param p : number of entries in the u_j array.
 * @param q_ij : matrix $q_{ij}$ wich results from $z_{ij}*u_{ij}$.
 * @param C_j : the sum of rows of $z_{ij}*u_{ij}$.
 * @param X : the datapoints matrix.
 * @param mu_j : the j_th cluster center.
 * @param M : result matrix @p M.
 * 
 * @code
 * template<typename T>
 * void
 * Kernel<T>::SMX(int n, int p,
 *                const T * q_ij,
 *                const T C_j,
 *                const int j,
 *                const T * X,
 *                const T * mu_j,
 *                T * M)
 * {
 *     int block_size = 512;
 *     int gridDim_y = (n+block_size-1)/block_size;
 * 
 *     dim3 num_blocks(p, gridDim_y);
 * 
 * 
 *     double time_sum = 0;
 *     int n_operations = 4;
 *     int n_memory_access = 5;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * @endcode
 * 
 * printf("Kernel Running ");
 * 

 * 
 * 
 * @code
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         __SMX<<<num_blocks, block_size>>>(n, p,
 *                                           q_ij,
 *                                           T(1./C_j),
 *                                           j,
 *                                           X,
 *                                           mu_j,
 *                                           M);
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * @endcode
 * 
 * if((s%10) == 0) printf(". ");
 * 

 * 
 * 
 * @code
 *     }
 * @endcode
 * 
 * printf("\n");
 * 
 * @code
 *     time_sum /= N_REPS;
 *     dt << "SMX_Kernel:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*p*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*p*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 * }
 * 
 * @endcode
 * 
 * This has to be at the end of file, because all functions have to be declared and their bodies defined before the
 * class can be explictly instantiated by the compiler.
 * 
 * @code
 * template class Kernel<float>;
 * 
 * #ifndef USE_SINGLE_PRECISION
 * template class Kernel<double>;
 * #endif
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef CUDADriver_STEP_7_H
 * #define CUDADriver_STEP_7_H
 * 
 * #include <deal.II/base/parameter_handler.h>
 * 
 * #include <lac/blas++.h>
 * 
 * #ifdef USE_KMEANS
 * #include <step-6/cuda_driver_step-6.h>
 * #include <step-6/cuda_driver_step-6.hh>
 * #include <step-6/Kmeans.h>
 * #endif
 * 
 * #include <step-7/cuda_kernel_wrapper_step-7.cu.h>
 * 
 * namespace step7 {
 * 
 * @endcode
 * 
 * 
 * <a name="ClassShohamEM"></a> 
 * <h3>Class: ShohamEM</h3>
 * 
 * 
 * declares the variables and methods needed for the Shoham EM algorithm
 * 
 * @code
 * template <typename T, typename blas>
 * class ShohamEM : public Kernel<T> {
 * public:
 * public:
 * 
 *     typedef typename blas_pp<T, blas>::blas_wrapper_type  BW;
 *     typedef typename blas_pp<T, blas>::FullMatrixAccessor FullMatrixAccessor;
 *     typedef typename blas_pp<T, blas>::Matrix             Matrix;
 *     typedef typename blas_pp<T, blas>::SubMatrix          SubMatrix;
 *     typedef typename blas_pp<T, blas>::MatrixSubCol       MatrixSubCol;
 *     typedef typename blas_pp<T, blas>::MatrixSubRow       MatrixSubRow;
 * 
 *     typedef typename blas_pp<T, blas>::SubColVector SubColVector;
 *     typedef typename blas_pp<T, blas>::SubRowVector SubRowVector;
 * 
 *     typedef typename blas_pp<T, blas>::Vector             Vector;
 *     typedef          transpose<Matrix>                    tr;
 *     typedef          transpose<SubMatrix>                 trSub;
 * 
 * @endcode
 * 
 * --- algorithm data
 * 
 * @code
 *     Matrix mu;
 *     int g;
 *     int g_min;
 *     std::vector<T> pi;
 *     std::vector<Matrix> Sigma;
 *     T v;
 *     double L;
 *     double L_max;
 *     double N;
 *     int p;
 *     Matrix u_ij;
 *     Matrix z_ij;
 *     Matrix P_ij;
 * 
 *     Matrix data;
 *     Matrix mean;
 *     int n_data_points;
 *     Matrix delta_ij;
 *     std::vector<Matrix> Sigma_inv;
 *     T gamma;
 *     Vector Delta;
 *     std::vector<Matrix> Sigma_dec;
 *     Vector pi_d;
 *     Matrix M;
 *     int iter;
 * @endcode
 * 
 * ---
 * 

 * 
 * 
 * @code
 *     std::vector<T> pi_opt;
 *     Matrix mu_opt;
 *     std::vector<Matrix> Sigma_opt;
 *     T v_opt;
 *     Matrix u_ij_opt;
 *     Matrix z_ij_opt;
 *     Matrix P_ij_opt;
 * 
 *     int it;
 * 
 *     ShohamEM(FullMatrixAccessor &ddata, FullMatrixAccessor &mmeans);
 * 
 *     void dump();
 * 
 *     void run_Shoham();
 *     void initialize_Shoham();
 *     void e_step_Shoham();
 *     void m_step_Shoham();
 *     void purge_step_Shoham();
 * 
 *     ~ShohamEM();
 * 
 * };
 * 
 * @endcode
 * 
 * 
 * <a name="ClassKernelTest"></a> 
 * <h3>Class: KernelTest</h3>
 * 
 * 
 * declares the variables and methods needed for the Kernel Tests
 * 
 * @code
 * template <typename T, typename blas>
 * class KernelTest : public Kernel<T> {
 * public:
 *     typedef typename blas_pp<T, blas>::blas_wrapper_type  BW;
 *     typedef typename blas_pp<T, blas>::FullMatrixAccessor FullMatrixAccessor;
 *     typedef typename blas_pp<T, blas>::Matrix             Matrix;
 *     typedef typename blas_pp<T, blas>::SubMatrix          SubMatrix;
 *     typedef typename blas_pp<T, blas>::MatrixSubCol       MatrixSubCol;
 *     typedef typename blas_pp<T, blas>::MatrixSubRow       MatrixSubRow;
 *     typedef typename blas_pp<T, blas>::Vector             Vector;
 *     typedef          transpose<Matrix>              tr;
 *     typedef          transpose<SubMatrix>              trSub;
 * 
 * 
 *     KernelTest(int n_data_points, int dim, int n_clusters);
 * 
 *     int n_data_points, n_clusters, dim;
 * 
 *     void run();
 * 
 *     ~KernelTest();
 * };
 * 
 * } // namespace step7 END
 * 
 * #endif // CUDADriver_STEP_7_H
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef CUDA_DRIVER_STEP_7_HH
 * #define CUDA_DRIVER_STEP_7_HH
 * #include <cuda_driver_step-7.h>
 * 
 * #include <cuda_kernel_wrapper_step-7.cu.h>
 * 
 * #include <cuda_kernel_wrapper_step-1.cu.h>
 * 
 * #include <base/CUDATimer.h>
 * 
 * #include <deal.II/base/timer.h>
 * 
 * #include <cmath>
 * #include <limits>
 * 
 * #include <lac/cublas_wrapper.hh>
 * 
 * #include <boost/math/special_functions/digamma.hpp>
 * #include <boost/math/special_functions/erf.hpp>
 * #include <boost/math/special_functions/beta.hpp>
 * 
 * #include <deal.II/numerics/histogram.h>
 * 
 * #include <iostream>
 * #include <fstream>
 * 
 * #include <QTime>
 * 
 * #include <deal.II/base/convergence_table.h>
 * #include <deal.II/base/table_handler.h>
 * #include <deal.II/base/timer.h>
 * 
 * #include <fstream>
 * 
 * @endcode
 * 
 * 
 * <a name="StructRightMSolve"></a> 
 * <h3>Struct: RightMSolve</h3>
 * Solves the equation ?X : AX = B
 * A,X,B Matrices
 * 
 * @code
 * template <typename T, typename blas>
 * struct RightMSolve{
 * 
 *     typedef bw_types::SubMatrixView<T, blas> SubMatrix;
 * 
 *     const SubMatrix &l;
 *     SubMatrix &r;
 * 
 *     RightMSolve(const SubMatrix & _l, SubMatrix & _r):l(_l),r(_r){}
 * 
 * @endcode
 * 
 * 
 * <a name="Operator"></a> 
 * <h4>Operator: = </h4>
 * overwrite the operator = for Matrix-Matrix Multiplication
 * 
 * @code
 *     RightMSolve<T, blas> & operator = (const SubMatrix & rhs)
 * @endcode
 * 
 * Argumens are consistance with the cublas Documentation
 * 
 * @code
 *     {
 * @endcode
 * 
 * copy the rightside to the resultsmatrix,so cublas will be able to overwrite
 * the memory area containing the rightside
 * with the solution.
 * 
 * @code
 *         r = rhs.matrix();
 * 
 *         std::cout << "? X : AX = B " << std::endl;
 *         char side = 'L', uplo = 'U', transa = 'N', diag = 'N';
 *         int m = rhs.r_end() - rhs.r_begin();
 *         int n = rhs.c_end() - rhs.c_begin();
 * 
 *         T alpha = 1.0;
 * 
 *         const T * const A = l.val();
 *         int lda = l.leading_dim();
 * 
 *         T *  B = r.val();
 *         int ldb = r.leading_dim();
 *         blas::trsm(  side, uplo, transa, diag,
 *                      m, n, alpha,
 *                      A, lda, B, ldb     );
 * 
 *         return *this;
 *     }
 * 
 * };
 * 
 * @endcode
 * 
 * 
 * <a name="StructRightMSolveTr"></a> 
 * <h3>Struct: RightMSolveTr</h3>
 * Solves the equation for $X$: $A^TX = B$ where $A, X, B$ are matrixes
 * 
 * @code
 * template <typename T, typename blas>
 * struct RightMSolveTr{
 * 
 *     typedef bw_types::SubMatrixView<T, blas> SubMatrix;
 *     typedef          transpose<SubMatrix>           tr;
 * 
 *     const SubMatrix &l;
 *     SubMatrix &r;
 * 
 *     RightMSolveTr(const tr & _l, SubMatrix & _r):l(_l.A),r(_r){}
 * 
 * @endcode
 * 
 * 
 * <a name="Operator"></a> 
 * <h4>Operator: = </h4>
 * overwrite the operator = for Matrix-Matrix Multiplication
 * 
 * @code
 *     RightMSolveTr<T, blas> & operator = (const SubMatrix & rhs)
 * @endcode
 * 
 * Argumens are consistance with the cublas Documentation
 * 

 * 
 * 
 * @code
 *     {
 * 
 * @endcode
 * 
 * copy the rightside to the resultsmatrix,so cublas will be able to overwrite
 * the memory area containing the rightside
 * with the solution.
 * 
 * @code
 *         r = rhs.matrix();
 * 
 *         std::cout << "? X : A^TX = B " << std::endl;
 *         char side = 'L', uplo = 'U', transa = 'T', diag = 'N';
 *         int m = rhs.r_end() - rhs.r_begin();
 *         int n = rhs.c_end() - rhs.c_begin();
 * 
 *         T alpha = 1.0;
 * 
 *         const T * const A = l.val();
 *         int lda = l.leading_dim();
 * 
 *         T * B = r.val();
 *         int ldb = r.leading_dim();
 *         cublas::trsm(  side, uplo, transa, diag,
 *                        m, n, alpha,
 *                        A, lda, B, ldb     );
 * 
 *         return *this;
 *     }
 * 
 * };
 * 
 * @endcode
 * 
 * 
 * <a name="Operator"></a> 
 * <h4>Operator: *</h4>
 * Overwrite the Operators * for Matrix-Matrix Multiplication
 * 
 * @code
 * template <typename T, typename blas>
 * inline RightMSolve<T, blas>  operator * (const bw_types::SubMatrixView<T, blas> &l, bw_types::SubMatrixView<T, blas> & r){
 *     return RightMSolve<T, blas>(l,r);
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Operator"></a> 
 * <h4>Operator: *</h4>
 * Overwrite the Operators * for Matrix-Matrix Multiplication
 * 
 * @code
 * template <typename T, typename blas>
 * inline RightMSolveTr<T, blas>  operator * (const typename RightMSolveTr<T, blas>::tr  &l,
 *                                            bw_types::SubMatrixView<T, blas> & r){
 *     return RightMSolveTr<T, blas>(l,r);
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="ConstructorShohamEM"></a> 
 * <h4>Constructor: ShohamEM</h4>
 * 
 * 
 * Creates an instance of the ShohamEm algorithm and initilize the data and mean matrices.
 * @param ddata : matrix containing the datapoints.
 * @param mmeans : matrix containing the starting means.
 * 
 * @code
 * template<typename T, typename blas>
 * step7::ShohamEM<T, blas>::ShohamEM(FullMatrixAccessor &ddata, FullMatrixAccessor &mmeans)
 *     :
 *       n_data_points(ddata.n_rows()), p(ddata.n_cols()),
 *       g(mmeans.n_rows()), data(ddata.n_rows(), ddata.n_cols())
 * {
 *     BW::Init();
 *     data = ddata;
 *     mean = mmeans;
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functiondump"></a> 
 * <h4>Function: dump</h4>
 * 
 * 
 * current data dump
 * 
 * @code
 * template<typename T, typename blas>
 * void step7::ShohamEM<T,blas>::dump()
 * {
 * 
 * @endcode
 * 
 * ---- save em means
 * 

 * 
 * 
 * @code
 *     std::fstream em;
 *     em.open("output/em.dat", std::ios::out);
 * 
 *     for(int i = 0; i < mu.n_rows(); i++){
 * 
 *         for(int j = 0; j < mu.n_cols(); j++)
 *             em << mu(i, j) << " ";
 * 
 *         em << "\n";
 * 
 *     }
 * 
 *     em.close();
 * 
 * @endcode
 * 
 * save Pi
 * 

 * 
 * 
 * @code
 *     std::fstream piout;
 *     piout.open("output/pi.dat", std::ios::out);
 * 
 *     for(int i = 0; i < pi.size(); i++)
 *         piout << pi[i] << "\n";
 * 
 *     piout.close();
 * 
 * @endcode
 * 
 * calculate ellipse parameters for each class
 * 

 * 
 * 
 * @code
 *     std::fstream ellipse;
 *     ellipse.open("output/ellipses.dat", std::ios::out);
 * 
 *     for(int k = 0; k < g; k++){
 * 
 *         double d;
 *         double eig[2], eigv[2];
 *         int perm[2];
 * 
 *         d = (Sigma[k](0, 0) - Sigma[k](1, 1)) / 2;
 * 
 *         eig[0] = (Sigma[k](0, 0) + Sigma[k](1, 1)) / 2 + sqrt(d * d + Sigma[k](0, 1) * Sigma[k](1, 0));
 *         eig[1] = (Sigma[k](0, 0) + Sigma[k](1, 1)) / 2 - sqrt(d * d + Sigma[k](0, 1) * Sigma[k](1, 0));
 * 
 *         if(eig[0] > eig[1] || true){
 *             perm[0] = 0;
 *             perm[1] = 1;
 *         }
 *         else{
 *             perm[0] = 1;
 *             perm[1] = 0;
 *         }
 * 
 *         std::cout << "\nEig: " << eig[perm[0]] << "\n";
 * 
 *         eigv[0] = eig[perm[0]] - Sigma[k](1, 1);
 *         eigv[1] = Sigma[k](1, 0);
 * 
 *         std::cout << "Eigv: " << eigv[0] / sqrt(eigv[0] * eigv[0] + eigv[1] * eigv[1]) << "..." << perm[0] << "\n";
 * 
 *         ellipse << mu(k, 0) << " " << mu(k, 1) << " " << / *sqrt(* /sqrt(sqrt(eig[perm[0]]))/ *)* / << " " << sqrt(sqrt(sqrt(eig[perm[1]])))
 *                 << " " <<  acos(eigv[0] / sqrt(eigv[0] * eigv[0] + eigv[1] * eigv[1])) * 180 / M_PI << "\n";
 * 
 *     }
 * 
 *     ellipse.close();
 * 
 * @endcode
 * 
 * assemble and execute gnuplot script
 * 

 * 
 * 
 * @code
 *     QString plot_waves, plot_pca, plot_ellipses, plot_hist, plot_process;
 * 
 *     FILE *gp = popen("/usr/local/lib/FEM/gnuplot/bin/gnuplot", "w");//open a pipe to gnuplot
 * 
 *     if(iter == 0){
 *         plot_waves += "!rm plot/ *\n";
 *     }
 * 
 *     plot_waves += "set terminal postscript landscape enhanced color solid linewidth 1.0 'Helvetica' 15\n"
 *             "set output 'plot/out" + QString::number(iter) + ".ps'\n"
 *             "set multiplot\n"
 *             "set xrange [0:128]\n"
 *             "set yrange [-1.5:1.5]\n"
 *             "set size 1,0.3\n"
 *             "set origin 0,0.8\n"
 *             "set key off\n";
 * 
 *     plot_pca += "reset\n"
 *             "set key off\n"
 *             "set size 1,0.7\n"
 *             "set origin 0,0.15\n"
 *             "set style data points\n"
 *             "set yrange []\n"
 *             "set xrange []\n";
 * 
 *     plot_hist = "set origin 0,-0.1\n"
 *             "set size 1,0.3\n"
 *             "set boxwidth 0.8\n"
 *             "set style fill solid\n"
 *             "set yrange [0:0.3]\n"
 *             "set xrange [-0.5:0.5]\n";
 * 
 *     for(int i = 1; i <= n_data_points; i+=1){
 * 
 *         int max = 0;
 * 
 *         for(int j = 1; j < z_ij.n_cols(); j++)
 *             if(z_ij(i-1, j) > z_ij(i-1, max))
 *                 max = j;
 * 
 *         int outlier = 1;
 *         if(u_ij(i-1, max) < 0.75) outlier = 4;
 *         if(u_ij(i-1, max) < 0.4) outlier = 6;
 * 
 *         outlier = 1;
 * 
 *         if(i == 1){
 *             plot_waves += "plot";
 *             plot_pca += "plot";
 *         }
 *         else{
 *             plot_waves += ",";
 *             plot_pca += ",";
 *         }
 * 
 *         plot_waves += " 'output/gnuplot_inputwaves.txt' u " + QString::number(i) +
 *                 "with lines lt 1 lc " + QString::number(max);
 * 
 *         plot_pca += " 'output/pca_out.txt' every ::" + QString::number(i - 1) + "::" + QString::number(i - 1) +
 *                 " lt " + QString::number(outlier) + " lc " + QString::number(max);
 * 
 *     }
 * 
 *     for(int i = 0; i < g; i++){
 * 
 *         plot_ellipses += ", 'output/ellipses.dat' every ::" + QString::number(i) + "::" + QString::number(i) +
 *                 " lc " + QString::number(i) + " with ellipses";
 * 
 *         if(i == 0)
 *             plot_hist += "plot";
 *         else
 *             plot_hist += ",";
 * 
 *         plot_hist += " 'output/pi.dat' every ::" + QString::number(i) + "::" + QString::number(i) +
 *                 " lc " + QString::number(i) + " with histograms";
 * 
 *     }
 * 
 *     plot_waves += "\n";
 *     fprintf(gp, plot_waves.toStdString().c_str());
 *     fflush(gp);
 * 
 *     fprintf(gp, plot_pca.toStdString().c_str());
 *     fflush(gp);
 * 
 *     fprintf(gp, ", 'output/em.dat' lt 3");
 * 
 *     plot_ellipses += "\n";
 *     fprintf(gp, plot_ellipses.toStdString().c_str());
 *     fflush(gp);
 * 
 *     plot_hist += "\n";
 *     fprintf(gp, plot_hist.toStdString().c_str());
 *     fflush(gp);
 * 
 *     plot_process = "!ps2pdf plot/out" + QString::number(iter) + ".ps plot/out" + QString::number(iter) + ".pdf\n"
 *             "!rm plot/out" + QString::number(iter) + ".ps\n"
 *             "!convert -density 300 plot/out" + QString::number(iter) + ".pdf plot/out" + QString::number(iter) + ".jpeg\n";
 * 
 *     int frames_per_image = 30; // for video only (set it to 1 if the video is not needed)
 * 
 *     for(int l = 1; l < frames_per_image; l++)
 *         plot_process += "!cp plot/out" + QString::number(iter) + ".jpeg plot/out" + QString::number(iter + l) + ".jpeg\n";
 * 
 *     fprintf(gp, plot_process.toStdString().c_str());
 *     fflush(gp);
 * 
 *     fprintf(gp, "unset multiplot\n");
 * 
 *     pclose(gp);
 * 
 *     iter += frames_per_image;
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionrun_Shoham"></a> 
 * <h4>Function: run_Shoham</h4>
 * 
 * 
 * runs the ShohamEM algorithm on given datapoints and means (new Version)
 * 
 * @code
 * template<typename T, typename blas>
 * void step7::ShohamEM<T, blas>::run_Shoham(){
 * 
 *     std::fstream file_like;
 *     file_like.open("output/likelihood.dat", std::ios::out);
 * 
 *     initialize_Shoham();
 * 
 *     dump();
 * 
 *     int iter = 0;
 *     T v_old;
 *     T L_old;
 * 
 *     while(g >= g_min){
 * 
 *         do{
 * 
 *             for(int o = 0; o < g; o++)
 *                 std::cout << pi[o];
 * 
 *             v_old = v;
 *             L_old = L;
 * 
 *             e_step_Shoham();
 * 
 *             m_step_Shoham();
 * 
 * @endcode
 * 
 * update P_ij
 * 
 * @code
 *             assemble_P_ij(gamma, v, p, n_data_points, g,
 *                           Delta.array().val(),
 *                           delta_ij.array().val(),
 *                           P_ij.array().val());
 * 
 * @endcode
 * 
 * update L
 * 
 * @code
 *             pi_d  = pi;
 * 
 *             Vector c_d;
 * 
 *             c_d = P_ij * pi_d;
 * 
 *             double log_1 = 0.;
 *             for(int j = 0; j< g; j++)
 *                 log_1 += std::log(n_data_points*pi[j]/12);
 * 
 *             double loglike_3 = (N/2)*log_1+(g/2)+(g*(N+1)/2);
 * 
 *             T sum;
 *             sum = sum_array(n_data_points, c_d.array().val(), true);
 * 
 *             L = sum - loglike_3;
 * 
 *             printf("###################loglike = %f ########################\n", L);
 * 
 *         } while(std::fabs(L_old - L) >= 0.1 || std::fabs(v - v_old) >= 0.01);
 * 
 *         dump();
 * 
 *         file_like << L << " "<<L_max <<" "<<v<< " " << std::fabs(L_old - L) << " "<<std::fabs(v - v_old)<<" "<<g<<"\n";
 *         flush(file_like);
 * 
 *         if(L >= L_max + L_max * 0.09 && g > g_min){
 * 
 *             L_max = L;
 * 
 *             pi_opt.resize(pi.size());
 *             std::copy(pi.begin(), pi.end(), pi_opt.begin());
 *             mu_opt = mu;
 *             Sigma_opt.resize(Sigma.size());
 *             std::copy(Sigma.begin(), Sigma.end(), Sigma_opt.begin());
 *             v_opt = v;
 *             u_ij_opt = u_ij;
 *             z_ij_opt = z_ij;
 *             P_ij_opt = P_ij;
 * 
 *             int min = 0;
 * 
 *             for(int i = 1; i < pi.size(); i++)
 *                 if(pi[i] < pi[min])
 *                     min = i;
 * 
 *             pi[min] = 0;
 * 
 *             purge_step_Shoham();
 * 
 *         }
 *         else
 *             break;
 * 
 *         iter++;
 * 
 *     }
 * 
 *     g++;
 *     pi = pi_opt;
 *     mu = mu_opt;
 *     Sigma = Sigma_opt;
 *     v = v_opt;
 *     u_ij = u_ij_opt;
 *     z_ij = z_ij_opt;
 *     P_ij = P_ij_opt;
 * 
 *     dump();
 * 
 *     file_like.close();
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functioninitialize_Shoham"></a> 
 * <h4>Function: initialize_Shoham</h4>
 * 
 * 
 * 
 * 
 * 
 * @code
 * template<typename T, typename blas>
 * void step7::ShohamEM<T, blas>::initialize_Shoham()
 * {
 *     it = 0;
 *     dealii::IdentityMatrix I(p);
 * 
 *     mu = mean; // cluster centers previous calculated by k-means
 * 
 * @endcode
 * 
 * g = g_max already set in constructor
 * 

 * 
 * 
 * @code
 *     g_min = 1;
 * 
 *     pi.resize(g);
 *     for(int i = 0; i < g; i++)
 *         pi[i] = 1./g;
 * 
 *     Sigma.resize(g);
 *     for(int i = 0; i < g; i++)
 *         Sigma[i] = I;
 * 
 *     v = 50;
 * 
 *     L_max = -1. * std::numeric_limits<T>::max(); // -infinity
 * 
 *     N = 0.1;//p;
 * 
 *     u_ij.reinit(n_data_points, g);
 *     z_ij.reinit(n_data_points, g);
 *     P_ij.reinit(n_data_points, g);
 * 
 *     delta_ij.reinit(n_data_points, g);
 *     Delta.reinit(g);
 * 
 *     Sigma_inv.resize(g);
 *     for(int i = 0; i < g; i++)
 *         Sigma_inv[i] = I;
 * 
 *     Sigma_dec.resize(g);
 *     for(int i = 0; i < g; i++)
 *         Sigma_dec[i] = I;
 * 
 *     M.reinit(n_data_points, p);
 * 
 *     iter = 0;
 * 
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functione_step_Shoham"></a> 
 * <h4>Function: e_step_Shoham</h4>
 * 
 * 
 * E Step of ShohamEM algorithm (update z and u)
 * 
 * @code
 * template<typename T, typename blas>
 * void step7::ShohamEM<T, blas>::e_step_Shoham()
 * {
 * 
 *     Matrix tmp(p, n_data_points);
 *     Vector c_d;
 *     std::vector<T> Delta_h(g);
 * 
 *     gamma = tgamma((v+p)/2) / (tgamma(v/2)*(std::pow((M_PI*v),(p/2.))));
 * 
 * @endcode
 * 
 * --- calculate delta_ij
 * 

 * 
 * 
 * @code
 *     for(int j=0;j<g;j++)
 *         Delta_h[j] = 1.;
 * 
 *     for(int k = 0; k < g; k++){
 * 
 *         subtract_array_from_matrix(n_data_points, p, k, g,
 *                                    data.array().val(),
 *                                    mu.array().val(),
 *                                    M.array().val());
 * 
 *         tmp = Sigma_inv[k] * tr(M);
 * 
 *         MM_scalar(n_data_points, p,
 *                   M.array().val(),
 *                   tmp.array().val(),
 *                   delta_ij.array().val() + k * n_data_points);
 *     }
 * 
 * 
 *     FullMatrixAccessor delta_ij_h(n_data_points, g, true);
 *     delta_ij_h = delta_ij;
 *     dealii::Vector<T> delta_ij_array(n_data_points * g);
 * 
 * 
 *     for(int i = 0; i <n_data_points; i++){ // transpose
 *         for(int j = 0; j < g; j++){
 *             delta_ij_array(j*n_data_points+i) = delta_ij(i,j);
 *         }
 *     }
 * 
 *     typename dealii::Vector<T>::const_iterator max_delta = std::max_element(delta_ij_array.begin(), delta_ij_array.end());
 *     typename dealii::Vector<T>::const_iterator min_delta = std::min_element(delta_ij_array.begin(), delta_ij_array.end());
 * 
 *     T interval = *max_delta - *min_delta;
 * 
 *     for(int i = 0; i<n_data_points*g; i++){
 *         delta_ij_array(i) -= *min_delta;
 *     }
 *     delta_ij_array /= interval;
 * 
 * 
 *     const unsigned int n_intervals = 50;
 *     dealii::Histogram Histogram;
 *     Histogram.evaluate(delta_ij_array, n_intervals, Histogram.linear);
 * 
 *     std::string delta_gram = "output/delta_gram.dat";
 *     std::ofstream dg(delta_gram.c_str());
 *     Histogram.write_gnuplot(dg);
 * 
 *     Delta = Delta_h;
 * 
 *     printf("gamma = %f\n", gamma);
 * 
 *     assemble_P_ij(gamma, v, p, n_data_points, g,
 *                   Delta.array().val(),
 *                   delta_ij.array().val(),
 *                   P_ij.array().val());
 * 
 *     pi_d  = pi;
 * 
 *     c_d = P_ij * pi_d;
 * 
 * @endcode
 * 
 * --- update z_ij
 * 
 * @code
 *     if(it++==0)
 *         assemble_z_ij(n_data_points, g,
 *                       P_ij.array().val(),
 *                       c_d.array().val(),
 *                       pi_d.array().val(),
 *                       z_ij.array().val());
 * 
 * 
 * 
 * @endcode
 * 
 * --- update u_ij
 * 

 * 
 * 
 * @code
 *     assemble_u_ij(v, p, n_data_points, g,
 *                   delta_ij.array().val(),
 *                   u_ij.array().val());
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionm_step_Shoham"></a> 
 * <h4>Function: m_step_Shoham</h4>
 * 
 * 
 * M Step of ShohamEM algorithm
 * 
 * @code
 * template<typename T, typename blas>
 * void step7::ShohamEM<T, blas>::m_step_Shoham()
 * {
 * 
 *     double digamma;
 *     Matrix P_ij_tmp = P_ij;
 *     Matrix z_ij_tmp = z_ij;
 *     Matrix q_ij(n_data_points, g);
 *     dealii::IdentityMatrix I(p);
 * 
 * @endcode
 * 
 * calculate y
 * 
 * @code
 *     digamma = boost::math::digamma((p+v)/2.);
 * 
 *     double y = calculate_y(n_data_points,
 *                            g,
 *                            v,
 *                            digamma,
 *                            z_ij.array().val(),
 *                            delta_ij.array().val(),
 *                            u_ij.array().val());
 * 
 *     y /= -n_data_points;
 * 
 * @endcode
 * 
 * y musst be greater than 1, because of the calculation of the DOF parameter.
 * 
 * @code
 *     y = std::max(y, 2.);
 *     printf("y = %f \n", y);
 * 
 * @endcode
 * 
 * update pi[j]
 * 
 * @code
 *     z_ij_tmp = z_ij;
 *     T max, sum_pi = 2., tmp_sum[g];
 * 
 *     for(int j = 0; j<g; j++){
 *         tmp_sum[j] = 0.;
 *     }
 * 
 *     int rel[g];
 * 
 *     for(int i = 0; i < g; i++) rel[i] = 0;
 * 
 *     for(int i = 0; i < n_data_points; i++){
 * 
 *         int max = 0;
 * 
 *         for(int j = 1; j < g; j++)
 *             if(z_ij(i, j) > z_ij(i, max))
 *                 max = j;
 * 
 *         rel[max]++;
 * 
 *     }
 * 
 *     int s = 0;
 * 
 *     for(int i = 0; i < g; i++)
 *         s+=rel[i];
 * 
 *     for(int i = 0; i < g; i++)
 *         pi[i] = ((double) rel[i]+1) / (double) s;
 * 
 *     for(int i = 0; i < g; i++)
 *         std::cout << pi[i] << ", ";
 * 
 *     purge_step_Shoham();
 * 
 *     std::cout << "\n\n";
 * 
 *     for(int i = 0; i < g; i++)
 *         std::cout << pi[i] << ", ";
 * 
 *     fflush(stdout);
 * 
 *     for(int i = 0; i < g; i++)
 *         if(pi[i] < 0.001) exit(0);
 * 
 * @endcode
 * 
 * update mu[j]
 * calculate sum(1-n)z_ij*u_ij
 * 
 * @code
 *     T zxu[g];
 *     for(int i = 0; i<g; i++){
 *         zxu[i] = 0.;
 *     }
 * 
 *     sum_array2D(n_data_points, g,
 *                 z_ij.array().val(),
 *                 u_ij.array().val(),
 *                 zxu);
 * 
 * @endcode
 * 
 * calculate z_ij*u_ij = q_ij
 * 
 * @code
 *     MxM(n_data_points, g,
 *         z_ij.array().val(),
 *         u_ij.array().val(),
 *         q_ij.array().val());
 * 
 * @endcode
 * 
 * update mu[j]
 * 
 * @code
 *     mu = tr(q_ij) * data;
 * 
 * 
 *     for(int i = 0; i<g; i++){
 *         MatrixSubRow mu_sub(mu,i,0);
 *         mu_sub *= 1./zxu[i];
 *     }
 * 
 * 
 * @endcode
 * 
 * update Sigma[j]
 * 
 * @code
 *     M.reinit(n_data_points, p);
 * 
 *     for(int j = 0; j<g; j++){
 *         SMX(n_data_points,
 *             p,
 *             q_ij.array().val(),
 *             zxu[j],
 *             j,
 *             data.array().val(),
 *             mu.array().val(),
 *             M.array().val()
 *             );
 * 
 *         Sigma[j] = tr(M)*M;
 *     }
 * 
 * 
 * @endcode
 * 
 * calculate Delta[j] and Sigma_inv[j]
 * 
 * @code
 *     int Matrix_Size = ((p+step1::DEFAULT_TILE_SIZE-1)/step1::DEFAULT_TILE_SIZE)*step1::DEFAULT_TILE_SIZE;
 *     for(int k=0; k < g; k++){
 * 
 * 
 * @endcode
 * 
 * cholesky decomposition
 * Save the small  Matrix in a Matrix TILE_SIZE $\times$ TILE_SIZE to use the Cholesky decomposition from Step_1
 * 
 * @code
 *         Matrix cholMatrix;
 *         dealii::IdentityMatrix eye(Matrix_Size);
 *         cholMatrix = eye;
 * 
 *         SubMatrix Sub_cholMatrix(cholMatrix, 0, p, 0 , p);
 * 
 *         Sub_cholMatrix = Sigma[k];
 * 
 *         step1::Kernels<T> fac_backend;
 * 
 *         fac_backend.cholesky.blackbox(cholMatrix.array().val(), cholMatrix.n_cols(), cholMatrix.leading_dim);
 * 
 *         SubMatrix Sub_Sigma(Sigma_dec[k], 0,0);
 * 
 *         Sub_Sigma = cholMatrix;
 * 
 *         T det = 0.;
 * 
 *         det = Det(p, Matrix_Size, Sub_cholMatrix.array().val());
 *         Delta.add(k, det);
 * 
 *         Matrix Identity = I;
 *         Matrix X(p,p);
 * 
 *         Sigma_inv[k] = Sigma_dec[k];
 *         const SubMatrix l_d(Sigma_dec[k],0,0);
 *         const SubMatrix rhs(Identity, 0,0);
 *         SubMatrix x_d(X,0,0);
 *         SubMatrix inv(Sigma_inv[k],0,0);
 * 
 *         RightMSolveTr<T, blas> L_X(trSub(l_d),x_d);
 *         L_X = rhs;
 * 
 *         RightMSolve<T, blas> L_B(l_d, inv);
 *         L_B = x_d;
 * 
 *     }
 * 
 * @endcode
 * 
 * update DOF parameter
 * 
 * @code
 *     double errfunc = boost::math::erf(0.6594*log(2.1971/(y+log(y)-1)));
 *     v = (2/(y+log(y)-1))+0.0416*(1+errfunc);
 * 
 *     printf("n_dofs_t_distro = %f\n", v);
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionpurge_step_Shoham"></a> 
 * <h4>Function: purge_step_Shoham</h4>
 * 
 * 
 * purges components in ShohamEM algorithm
 * 
 * @code
 * template<typename T, typename blas>
 * void step7::ShohamEM<T, blas>::purge_step_Shoham()
 * {
 * 
 *     for(int i = 0; i < pi.size(); i++){
 * 
 *         if(pi[i] < 0.0001 && g>g_min){
 * 
 * @endcode
 * 
 * update pi
 * 
 * @code
 *             pi[i] = pi[g-1];
 *             pi.resize(g-1);
 * 
 * @endcode
 * 
 * update sigma
 * 
 * @code
 *             Sigma[i] = Sigma[g-1];
 *             Sigma.resize(g-1);
 * 
 * @endcode
 * 
 * update P_ij
 * 
 * @code
 *             MatrixSubCol subPij_target(P_ij, 0, i);
 *             MatrixSubCol subPij_source(P_ij, 0, g-1);
 *             subPij_target = subPij_source;
 * 
 *             Matrix new_Pij(n_data_points, g-1);
 * 
 *             for(int j = 0; j < g-1; j++){
 * 
 *                 MatrixSubCol to(new_Pij, 0, j);
 *                 MatrixSubCol from(P_ij, 0, j);
 * 
 *                 to = from;
 * 
 *             }
 * 
 *             P_ij = new_Pij;
 * 
 * @endcode
 * 
 * update delta_ij
 * 
 * @code
 *             MatrixSubCol subDeltaij_target(P_ij, 0, i);
 *             MatrixSubCol subDeltaij_source(P_ij, 0, g-1);
 *             subDeltaij_target = subDeltaij_source;
 * 
 *             Matrix new_Deltaij(n_data_points, g-1);
 * 
 *             for(int j = 0; j < g-1; j++){
 * 
 *                 MatrixSubCol to(new_Deltaij, 0, j);
 *                 MatrixSubCol from(delta_ij, 0, j);
 * 
 *                 to = from;
 * 
 *             }
 * 
 *             delta_ij = new_Deltaij;
 * 
 * @endcode
 * 
 * update mu;
 * 
 * @code
 *             MatrixSubRow subMu_target(mu, i, 0);
 *             MatrixSubRow subMu_source(mu, g-1, 0);
 *             subMu_target = subMu_source;
 * 
 *             Matrix new_mu(g-1, p);
 * 
 *             for(int j = 0; j < g-1; j++){
 * 
 *                 MatrixSubRow to(new_mu, j, 0);
 *                 MatrixSubRow from(mu, j, 0);
 * 
 *                 to = from;
 * 
 *             }
 * 
 *             mu = new_mu;
 * 
 * @endcode
 * 
 * update u_ij
 * 
 * @code
 *             MatrixSubCol subUij_target(u_ij, 0, i);
 *             MatrixSubCol subUij_source(u_ij, 0, g-1);
 *             subUij_target = subUij_source;
 * 
 *             Matrix new_Uij(n_data_points, g-1);
 * 
 *             for(int j = 0; j < g-1; j++){
 * 
 *                 MatrixSubCol to(new_Uij, 0, j);
 *                 MatrixSubCol from(u_ij, 0, j);
 * 
 *                 to = from;
 * 
 *             }
 * 
 *             u_ij = new_Uij;
 * 
 * @endcode
 * 
 * update z_ij
 * 
 * @code
 *             MatrixSubCol subZij_target(z_ij, 0, i);
 *             MatrixSubCol subZij_source(z_ij, 0, g-1);
 *             subZij_target = subZij_source;
 * 
 *             Matrix new_Zij(n_data_points, g-1);
 * 
 *             for(int j = 0; j < g-1; j++){
 * 
 *                 MatrixSubCol to(new_Zij, 0, j);
 *                 MatrixSubCol from(z_ij, 0, j);
 * 
 *                 to = from;
 * 
 *             }
 * 
 *             z_ij = new_Zij;
 * 
 * @endcode
 * 
 * update Delta
 * 
 * @code
 *             Delta.set(i, Delta(g-1));
 *             Vector newDelta;
 *             newDelta.reinit(g-1);
 * 
 *             for(int j = 0; j < g-1; j++)
 *                 newDelta.set(j, Delta(j));
 * 
 *             Delta = newDelta;
 * 
 * @endcode
 * 
 * update Sigma_inv
 * 
 * @code
 *             Sigma_inv[i] = Sigma_inv[g-1];
 *             Sigma_inv.resize(g-1);
 * 
 * @endcode
 * 
 * update Sigma_dec
 * 
 * @code
 *             Sigma_dec[i] = Sigma_dec[g-1];
 *             Sigma_dec.resize(g-1);
 * 
 *             g--;
 * 
 *             if(pi[i] == 0) // copied one is zero
 *                 i--;
 * 
 *             purge_step_Shoham();
 * 
 *             break;
 * 
 *         }
 * 
 *     }
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="DestructorShohamEM"></a> 
 * <h4>Destructor : ~ShohamEM</h4>
 * 
 * 
 * Free the memory taken by ShohamEM and shuts down the blas library.
 * 
 * @code
 * template<typename T, typename blas>
 * step7::ShohamEM<T, blas>::~ShohamEM()
 * {
 *     BW::Shutdown();
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="ConstructorKernelTest"></a> 
 * <h4>Constructor: KernelTest</h4>
 * 
 * 
 * Initilize the Kernel performance tests
 * @param n_data_points : number of datapoints.
 * @param dim : dimension (number of columns).
 * @param n_clusters : number of clusters.
 * 
 * @code
 * template<typename T, typename blas>
 * step7::KernelTest<T, blas>::KernelTest(int n_data_points, int dim, int n_clusters)
 *     :
 *       n_data_points(n_data_points), dim(dim),
 *       n_clusters(n_clusters)
 * {
 *     BW::Init();
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionrun"></a> 
 * <h4>Function: run</h4>
 * 
 * 
 * runs the KernelTest
 * 
 * @code
 * template<typename T, typename blas>
 * void step7::KernelTest<T, blas>::run()
 * {
 *     int Matrix_Size = ((dim+step1::DEFAULT_TILE_SIZE-1)/step1::DEFAULT_TILE_SIZE)*step1::DEFAULT_TILE_SIZE;
 * 
 *     FullMatrixAccessor Test_Matrix(n_data_points, dim, true);
 *     FullMatrixAccessor Test_Matrix2(n_clusters, dim, true);
 *     FullMatrixAccessor Test_Matrix3(n_data_points, dim, true);
 *     FullMatrixAccessor Test_Matrix4(dim, n_data_points, true);
 * 
 *     FullMatrixAccessor Test_Matrix_ij(n_data_points, n_clusters, true);
 *     FullMatrixAccessor Test_Matrix_ij2(n_data_points, n_clusters, true);
 *     FullMatrixAccessor Test_Matrix_ij3(n_data_points, n_clusters, true);
 * 
 *     std::vector<T> Delta_h(n_clusters);
 *     std::vector<T> c_h(n_data_points);
 *     std::vector<T> pi_h(n_clusters);
 *     Vector Delta, c, pi;
 * 
 *     T zxu[n_clusters];
 *     for(int i = 0; i<n_clusters; i++){
 *         zxu[i] = 0.;
 *     }
 * 
 *     for(int i = 0; i< n_data_points; i++){
 *         for(int j=0; j<dim; j++){
 *             Test_Matrix(i,j) = (i+1)/10.;
 *             Test_Matrix3(i,j) = 0.;
 *             Delta_h[j] = 1.;
 *             c_h[i] = 1.;
 *             pi_h[j] = 1/n_clusters;
 *         }
 *     }
 * 
 *     for(int i = 0; i< n_clusters; i++){
 *         for(int j=0; j<dim; j++){
 *             Test_Matrix2(i,j) = (i+1)/10.;
 *         }
 *     }
 * 
 *     for(int i = 0; i< dim; i++){
 *         for(int j=0; j<n_data_points; j++){
 *             Test_Matrix4(i,j) = (i+1)/10.;
 *         }
 *     }
 * 
 *     for(int i = 0; i< n_data_points; i++){
 *         for(int j=0; j<n_clusters; j++){
 *             Test_Matrix_ij(i,j) = (i+1)/10.;
 *             Test_Matrix_ij2(i,j) = (i+1)/10.;
 *             Test_Matrix_ij3(i,j) = (i+1)/10.;
 *         }
 *     }
 * 
 *     Matrix Test, Test2, Test3, Test4, Test_ij, Test_ij2, Test_ij3;
 *     Test = Test_Matrix;
 *     Test2 = Test_Matrix2;
 *     Test3 = Test_Matrix3;
 *     Test4 = Test_Matrix4;
 *     Test_ij = Test_Matrix_ij;
 *     Test_ij2 = Test_Matrix_ij2;
 *     Test_ij3 = Test_Matrix_ij3;
 *     Delta = Delta_h;
 *     c = c_h;
 *     pi = pi_h;
 * 
 *     subtract_array_from_matrix(n_data_points, dim, 1, n_clusters,
 *                                Test.array().val(),
 *                                Test2.array().val(),
 *                                Test3.array().val());
 * 
 *     printf("dim = %d\n", dim);
 * 
 * 
 *     MM_scalar(n_data_points, dim,
 *               Test3.array().val(),
 *               Test4.array().val(),
 *               Test.array().val());
 * 
 *     Det(dim, Matrix_Size, Test3.array().val());
 * 
 *     assemble_P_ij(0.5, 50, dim, n_data_points, n_clusters,
 *                   Delta.array().val(),
 *                   Test_ij.array().val(),
 *                   Test_ij2.array().val());
 * 
 *     assemble_u_ij(50, dim, n_data_points, n_clusters,
 *                   Test_ij.array().val(),
 *                   Test_ij2.array().val());
 * 
 *     assemble_z_ij(n_data_points, n_clusters,
 *                   Test_ij.array().val(),
 *                   c.array().val(),
 *                   pi.array().val(),
 *                   Test_ij2.array().val());
 * 
 *     calculate_y(n_data_points, n_clusters, 0.5, 13,
 *                 Test_ij.array().val(),
 *                 Test_ij2.array().val(),
 *                 Test_ij3.array().val());
 * 
 *     sum_array(n_data_points, c.array().val(), false);
 * 
 *     sum_array2D(n_data_points, n_clusters,
 *                 Test_ij.array().val(),
 *                 Test_ij2.array().val(),
 *                 zxu);
 * 
 *     MxM(n_data_points, n_clusters,
 *         Test_ij.array().val(),
 *         Test_ij2.array().val(),
 *         Test_ij3.array().val());
 * 
 *     SMX(n_data_points, dim,
 *         Test_ij.array().val(),
 *         zxu[2],
 *         2,
 *         Test.array().val(),
 *         Test2.array().val(),
 *         Test3.array().val());
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="DestructorKernelTest"></a> 
 * <h4>Destructor : ~KernelTest</h4>
 * 
 * 
 * Free the memory taken by KernelTest and shuts down the blas library.
 * 
 * @code
 * template<typename T, typename blas>
 * step7::KernelTest<T, blas>::~KernelTest()
 * {
 *     BW::Shutdown();
 * }
 * #endif // CUDA_DRIVER_STEP_7_HH
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * @endcode
 * 
 * STL header
 * 
 * @code
 * #include <iostream>
 * #include <vector>
 * #include <fstream>
 * 
 * #include <QString>
 * 
 * @endcode
 * 
 * Driver for the GPU-Part
 * 
 * @code
 * #include <cuda_driver_step-7.h>
 * 
 * @endcode
 * 
 * deal.II-Komponenten
 * 
 * @code
 * #include <deal.II/lac/matrix_out.h>
 * 
 * @endcode
 * 
 * cublas-Wrapper-Classes.<br>
 * attach all the needed header-files.
 * 
 * @code
 * #include <cuda_driver_step-7.hh>
 * 
 * #include <waveforms_generator.h>
 * 
 * @endcode
 * 
 * Reuse GlobalData and the Params classes from @ref step_4 "step-4".
 * 
 * @code
 * #include <step-4.hh>
 * 
 * @endcode
 * 
 * Use NIPALS from @ref step_5 "step-5"
 * 
 * @code
 * #include <cuda_driver_step-5.h>
 * #include <cuda_driver_step-5.hh>
 * 
 * namespace step7 {
 * 
 * @endcode
 * 
 * 
 * <a name="ClassSimParams"></a> 
 * <h3>Class: SimParams</h3>
 * contains the the parameter needed to run the simulation
 * 
 * @code
 * template<typename T>
 * struct SimParams : public step4::QRTestUIParams
 *         #ifdef USE_KMEANS
 *         , step6::SimParam
 *         #endif
 * 
 * {
 * 
 *     typedef step4::QRTestUIParams Base1;
 * 
 * #ifdef USE_KMEANS
 *     typedef step6::SimParam Base2;
 * #endif
 * 
 *     SimParams() : Base1()
 *   #ifdef USE_KMEANS
 *       , Base2()
 *   #endif
 *     {}
 * 
 *     static void declare(dealii::ParameterHandler & prm);
 * 
 *     void get(dealii::ParameterHandler & prm);
 * 
 * @endcode
 * 
 * waveforms parameter
 * 
 * @code
 *     dealii::FullMatrix<T> waveforms;
 *     dealii::FullMatrix<T> check_waves;
 * 
 *     int n_forms, n_rows, n_cols;
 *     const static int forms = 3;
 *     T noise;
 *     std::vector<double> waves_parameters;
 *     std::vector<int> n_waves;
 *     QString number_of_waves;
 *     std::string output_folder;
 * 
 *     unsigned int n_components, max_iter;
 * 
 *     double ev_tol;
 * 
 * public:
 * 
 *     void generate_waveforms();
 * 
 * };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ClassPCA"></a> 
 * <h3>Class: PCA</h3>
 * 
 * 
 * contains all the needed variables and function for the PCA
 * 
 * @code
 * template <typename T, typename blas> class PCA {
 * 
 * public:
 * 
 *     typedef typename blas_pp<T, blas>::FullMatrixAccessor FullMatrixAccessor;
 *     typedef typename blas_pp<T, blas>::Matrix       Matrix;
 *     typedef          transpose<Matrix>              tr;
 * 
 *     typedef typename blas_pp<T, blas>::SubMatrix    SubMatrix;
 *     typedef typename blas_pp<T, blas>::MatrixSubCol MatrixSubCol;
 *     typedef typename blas_pp<T, blas>::Vector       Vector;
 * 
 * 
 *     PCA(dealii::ParameterHandler & prm);
 *     void run();
 * 
 * 
 *     Matrix Q_x_scores;
 *     Matrix d_scores;
 *     Matrix d_loadings;
 * 
 * private:
 *     void factorize(dealii::FullMatrix<T> &A);
 * 
 *     void check_results(const dealii::FullMatrix<T> & A,
 *                        const step4::CudaQRDecomposition<T, blas>& QRf,
 *                        T elapsed_time);
 * 
 *     void save_results();
 * 
 *     dealii::FullMatrix<T> Q, B, P_t, H;
 * 
 *     SimParams<T> params;
 * 
 *     dealii::ConvergenceTable results_table;
 * 
 *     unsigned int n_successful_measurements;
 * };
 * 
 * @endcode
 * 
 * 
 * <a name="ClassMyFancySimulation"></a> 
 * <h3>Class: MyFancySimulation</h3>
 * 
 * 
 * This class controls the actual simulation for the physical problem.
 * 
 * @code
 * template<typename T>
 * class MyFancySimulation {
 * 
 * public:
 * 
 *     typedef typename blas_pp<T, blas>::FullMatrixAccessor FullMatrixAccessor;
 *     typedef typename blas_pp<T, blas>::Matrix       Matrix;
 *     typedef          transpose<Matrix>              tr;
 * 
 *     typedef typename blas_pp<T, blas>::SubMatrix    SubMatrix;
 *     typedef typename blas_pp<T, blas>::MatrixSubCol MatrixSubCol;
 *     typedef typename blas_pp<T, blas>::Vector       Vector;
 * 
 *     MyFancySimulation(std::string prm_filename);
 * 
 *     void run();
 * 
 *     void run_pca(dealii::ParameterHandler &prm_handler);
 * 
 *     void run_cluster_analysis();
 * 
 * private:
 * 
 *     SimParams<T> params;
 * 
 *     FullMatrixAccessor data;
 *     FullMatrixAccessor loadings;
 *     FullMatrixAccessor  scores;
 * };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ConstructorPCA"></a> 
 * <h4>Constructor: PCA</h4>
 * 
 * 
 * Init the PCA class the get the needed parameter from the .prm file
 * @param prm : dealii parameter handler.
 * 
 * @code
 * template <typename T, typename blas>
 * step7::PCA<T, blas>::PCA(dealii::ParameterHandler &prm)
 *     :
 *       n_successful_measurements(0)
 * {
 *     params.get(prm);
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionrun"></a> 
 * <h4>Function: run</h4>
 * 
 * 
 * generates Data for the PCA and calls the factorize function
 * 
 * @code
 * template <typename T, typename blas>
 * void step7::PCA<T, blas>::run() {
 * 
 *     params.generate_waveforms();
 * 
 * @endcode
 * 
 * write the calculated means and covariances to the files means.txt and sig.txt respectively
 * 
 * @code
 *     std::string wforms = "output/wforms.txt";
 *     std::ofstream wf(wforms.c_str());
 *     params.waveforms.print(wf, 12 / *width; avoids that two numbers appear as one in the output* /);
 * 
 *     factorize(params.waveforms);
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionfactorize"></a> 
 * <h4>Function: factorize</h4>
 * 
 * 
 * runs the QR-Decomposition on the Data-Matrix followes by the PCA using the nipals Algorithm
 * @param A : input Matrix.
 * 
 * @code
 * template <typename T, typename blas>
 * void
 * step7::PCA<T, blas>::factorize(dealii::FullMatrix<T> &A)
 * {
 * 
 *     std::cout  << std::endl << " ---------- FACTORIZE ----------" << std::endl;
 * 
 *     step4::CudaQRDecomposition<T, blas> QR;
 * 
 *     QR.householder(A);
 * 
 *     std::cout  << std::endl << " ---------- DONE ----------" << std::endl;
 * 
 *     FullMatrixAccessor R_tmp;
 *     R_tmp = QR.R();
 * 
 * 
 *     dealii::FullMatrix<T> R_tmp_h(R_tmp.n_rows(), R_tmp.n_cols());
 * 
 *     std::ofstream R_out("output/R_mea.dat");
 * 
 *     std::cout<< "R dimensions = " << R_tmp.n_rows() << " X " << R_tmp.n_cols() << std::endl;
 *     for (unsigned int r = 0; r < R_tmp_h.n_rows(); ++r)
 *     {
 *         for(unsigned int c = 0; c < R_tmp_h.n_cols(); ++c)
 *         {
 *             R_tmp_h(r, c) = R_tmp(r, c);
 *             R_out << R_tmp_h(r, c) << "\t";
 *         }
 *         R_out  << std::endl;
 *     }
 *     std::cout  << std::endl;
 * 
 *     std::cout  << std::endl << " ---------- NIPALS ----------" << std::endl;
 * 
 *     step5::CudaNIPALS<T, blas> nipals(params.n_components,
 *                                       -1/ *params.max_iter* /,
 *                                       params.ev_tol);
 * 
 *     nipals.get_pca(R_tmp_h);
 * 
 *     std::cout  << std::endl << " ---------- DONE ----------" << std::endl;
 * 
 *     nipals.lambda.print(std::cout);
 * 
 * 
 *     Q_x_scores.reinit(params.n_rows,params.n_components);
 * 
 *     SubMatrix sub_Q(const_cast<Matrix &>(QR.Q()),0,params.n_rows,0,params.n_cols);
 * 
 *     SubMatrix scores(const_cast<Matrix &>(nipals.scores()),0,0);
 * 
 *     SubMatrix Q_x_s(const_cast<Matrix &>(Q_x_scores),0,0);
 * 
 * 
 *     d_scores.reinit(nipals.scores().n_rows(), nipals.scores().n_cols());
 *     d_scores = nipals.scores();
 * 
 *     d_loadings.reinit(nipals.loads().n_rows(),nipals.loads().n_cols());
 *     d_loadings = nipals.loads();
 * 
 *     Matrix Q_tmp, S_tmp, P_tmp;
 *     Q_tmp = QR.Q();
 *     S_tmp = nipals.scores();
 *     P_tmp = Q_tmp * S_tmp;
 * 
 *     Q_x_scores = P_tmp;
 * 
 *     std::ofstream PCA_out("output/pca_out.txt");
 *     for(unsigned int r = 0; r < P_tmp.n_rows(); r++){
 *         for(unsigned int c = 0; c < P_tmp.n_cols(); c++)
 *             PCA_out<<P_tmp(r,c)<<"\t";
 *         PCA_out<<"\n";
 *     }
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functiondeclare"></a> 
 * <h4>Function: declare</h4>
 * 
 * 
 * Declares the input parameter for the dealII Parameter Handler
 * @param prm : dealii parameter handler.
 * 
 * @code
 * template<typename T>
 * void step7::SimParams<T>::declare(dealii::ParameterHandler &prm)
 * {
 *     Base1::declare(prm);
 * #ifdef USE_KMEANS
 *     Base2::declare(prm);
 * #endif
 * 
 *     prm.enter_subsection("Waveforms parameters");
 * 
 *     prm.declare_entry("Number of waves", "300,300,300",
 *                       dealii::Patterns::Anything(),
 *                       "Number of Waves");
 * 
 *     prm.declare_entry("Noise intensity", "0.2",
 *                       dealii::Patterns::Double(0.,0.9),
 *                       "Noise intensity");
 * 
 *     prm.declare_entry("Number of columns", "128",
 *                       dealii::Patterns::Integer(),
 *                       "Number of columns");
 * 
 *     prm.declare_entry("Center of Wave 1", "0.25",
 *                       dealii::Patterns::Double(0.,0.9),
 *                       "Center of Wave 1");
 * 
 *     prm.declare_entry("Sigma of Wave 1", "0.09",
 *                       dealii::Patterns::Double(0.01,0.09),
 *                       "Sigma of Wave 1 (standard deviation of gaussian)");
 * 
 *     prm.declare_entry("Center of Wave 2", "0.85",
 *                       dealii::Patterns::Double(0., 1.),
 *                       "Center of Wave 2");
 * 
 *     prm.declare_entry("Center 1 of Wave 3", "0.6",
 *                       dealii::Patterns::Double(0.,0.9),
 *                       "Center 1 of Wave 3");
 * 
 *     prm.declare_entry("Center 2 of Wave 3", "0.65",
 *                       dealii::Patterns::Double(0.,0.9),
 *                       "Center 2 of Wave 3");
 * 
 *     prm.declare_entry("Sigma 1 of Wave 3", "0.05",
 *                       dealii::Patterns::Double(0.01,0.09),
 *                       "Sigma 1 of Wave 3 (Lorenz)");
 * 
 *     prm.declare_entry("Sigma 2 of Wave 3", "0.06",
 *                       dealii::Patterns::Double(0.01,0.09),
 *                       "Sigma 2 of Wave 3");
 * 
 *     prm.declare_entry ("Output folder",
 *                        "output",
 *                        dealii::Patterns::Anything(),
 *                        "output folder location");
 * 
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("NIPALS parameters");
 * 
 *     prm.declare_entry("Number of components", "2",
 *                       dealii::Patterns::Integer(),
 *                       "Number of components");
 * 
 *     prm.declare_entry("Maximal number of iterations", "1000",
 *                       dealii::Patterns::Integer(),
 *                       "Maximal number of iterations. "
 *                       "If -1 is given termination of NIPALS iterations "
 *                       "solely depends on the given tolerance");
 * 
 *     prm.declare_entry("Tolerance for eigenvalues",
 *                       "1e-5",
 *                       dealii::Patterns::Double(),
 *                       "Once the relative change of an eigenvalue in the inner "
 *                       "power iteration is les than the given value iteration "
 *                       "stops and the next principal component is computed.");
 * 
 *     prm.leave_subsection();
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionget"></a> 
 * <h4>Function: get</h4>
 * 
 * 
 * Reads the input parameters from the .prm file
 * @param prm : dealii parameter handler.
 * 
 * @code
 * template<typename T>
 * void step7::SimParams<T>::get(dealii::ParameterHandler &prm)
 * {
 * 
 *     this->Base1::get(prm);
 * 
 * #ifdef USE_KMEANS
 *     this->Base2::get(prm);
 * #endif
 * 
 *     prm.enter_subsection("Waveforms parameters");
 * 
 *     number_of_waves = QString(prm.get("Number of waves").c_str());
 * 
 *     noise = prm.get_double("Noise intensity");
 * 
 *     n_cols = prm.get_integer("Number of columns");
 * 
 *     waves_parameters.push_back(prm.get_double("Center of Wave 1"));
 * 
 *     waves_parameters.push_back(prm.get_double("Sigma of Wave 1"));
 * 
 *     waves_parameters.push_back(prm.get_double("Center of Wave 2"));
 * 
 *     waves_parameters.push_back(prm.get_double("Center 1 of Wave 3"));
 * 
 *     waves_parameters.push_back(prm.get_double("Center 2 of Wave 3"));
 * 
 *     waves_parameters.push_back(prm.get_double("Sigma 1 of Wave 3"));
 * 
 *     waves_parameters.push_back(prm.get_double("Sigma 2 of Wave 3"));
 * 
 *     output_folder = prm.get ("Output folder");
 * 
 *     QDir dir("plot");
 * 
 *     if(!dir.exists()){
 *         bool mkdir = dir.mkpath("./");
 *         AssertThrow(mkdir, dealii::ExcMessage("creating a new Folder failed"));
 *     }
 * 
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("NIPALS parameters");
 * 
 *     this->n_components = prm.get_integer("Number of components");
 * 
 *     this->max_iter = prm.get_integer("Maximal number of iterations");
 * 
 *     this->ev_tol   = prm.get_double("Tolerance for eigenvalues");
 * 
 *     prm.leave_subsection();
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functiongenerate_waveforms"></a> 
 * <h4>Function: generate_waveforms</h4>
 * 
 * 
 * Generates a test matrix for the algorithm
 * The data points are read from the .prm file
 * 
 * @code
 * template<typename T>
 * void step7::SimParams<T>::generate_waveforms()
 * {
 *     QStringList number_list = number_of_waves.split(",", QString::SkipEmptyParts);
 * 
 *     if(number_list.size()<forms)
 *     {
 *         std::cout<<"more waves numbers needed!"<<std::endl;
 *         exit(1);
 *     }
 * 
 *     QStringList::iterator e = number_list.begin();
 *     QStringList::iterator end = number_list.end();
 * 
 *     int i=0, sum = 0;
 *     n_waves.resize(forms);
 *     n_forms = 0;
 * 
 * @endcode
 * 
 * Read the number of waves from the .prm file
 * using split to seperate the string into single values
 * 
 * @code
 *     for(; e!= end; ++e, ++i)
 *     {
 *         if(i==forms) break;
 *         n_waves[i] = e->toInt();
 *     }
 * 
 *     for(int i = 0; i<forms;i++){
 *         sum += n_waves[i];
 *         if(n_waves[i]!=0) n_forms++;
 *     }
 * 
 *     n_rows = sum;
 * 
 *     waveforms.reinit(n_rows,n_cols);
 *     check_waves.reinit(n_forms,n_cols);
 *     WaveForms::generate_waveforms(n_rows, n_cols, waveforms, n_waves, noise, waves_parameters);
 *     WaveForms::normalize_waveforms(n_rows, n_cols, waveforms);
 *     WaveForms::original_waveforms(n_cols,n_forms,check_waves, waves_parameters, n_waves);
 * 
 * @endcode
 * 
 * Write the waveforms to File
 * 
 * @code
 *     dealii::MatrixOut matrix_out;
 *     std::string wave_forms = output_folder+"/waveforms.txt";
 *     std::ofstream waf(wave_forms.c_str());
 *     matrix_out.build_patches(waveforms, "waveforms");
 *     matrix_out.write_gnuplot(waf);
 * 
 * @endcode
 * 
 * Write the originalforms to File
 * 
 * @code
 *     dealii::MatrixOut original_out;
 *     std::string original_forms = output_folder+"/originalforms.txt";
 *     std::ofstream oaf(original_forms.c_str());
 *     original_out.build_patches(check_waves, "originalforms");
 * @endcode
 * 
 * original_out.write_gnuplot(oaf);
 * 

 * 
 * 
 * @code
 * }
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="ConstructorMyFancySimulation"></a> 
 * <h4>Constructor: MyFancySimulation</h4>
 * 
 * 
 * The constructor initialize the simulation and calls the dealII parameter Handler to get the parameter from the .prm file
 * @param prm_filename : name of the parameter file.
 * 
 * @code
 * template<typename T>
 * step7::MyFancySimulation<T>::MyFancySimulation(std::string prm_filename)
 * {
 *     dealii::ParameterHandler param_handler;
 *     SimParams<T>::declare(param_handler);
 * 
 *     param_handler.read_input(prm_filename);
 *     params.get(param_handler);
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Funktionrun"></a> 
 * <h4>Funktion: run</h4>
 * 
 * 
 * The real simulation will be started with the run() function
 * 
 * @code
 * template<typename T>
 * void step7::MyFancySimulation<T>::run()
 * {
 * 
 * 
 * }
 * 
 * template<typename T>
 * void step7::MyFancySimulation<T>::run_pca(dealii::ParameterHandler &prm_handler)
 * {
 * 
 *     std::cout << "Householder-QR mit cublas :\n"
 *               << "------------------------------------------------" << std::endl;
 * 
 *     PCA<T, blas> cublas_pca(prm_handler);
 * 
 *     cublas_pca.run();
 *     data = cublas_pca.Q_x_scores;
 *     loadings = cublas_pca.d_loadings;
 *     scores = cublas_pca.d_scores;
 *     std::cout << "\nHouseholder-QR mit cublas DONE \n"
 *               << "------------------------------------------------\n" << std::endl;
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Funktionrun_cluster_analysis"></a> 
 * <h4>Funktion: run_cluster_analysis</h4>
 * 
 * 
 * runs the clusteranalysis using The Shoham Algorithm or K-Means
 * 
 * @code
 * template<typename T>
 * void step7::MyFancySimulation<T>::run_cluster_analysis()
 * {
 * 
 *     std::cout  << std::endl << " ---------- START SIMULATION ----------" << std::endl;
 * 
 *     std::cout  << std::endl << " ---------- GENERATE WAVEFORMS ----------" << std::endl;
 * 
 *     params.generate_waveforms();
 * 
 *     std::cout  << std::endl << " ---------- DONE ----------" << std::endl;
 * 
 *     std::cout  << std::endl;
 * 
 *     FullMatrixAccessor X(params.waveforms,true);
 * 
 *     FullMatrixAccessor Y(params.check_waves, true);
 * 
 *     dealii::FullMatrix<T> data_matrix2(params.n_rows, params.n_cols);
 *     dealii::FullMatrix<T> data_matrix2t(params.n_cols, params.n_rows);
 *     X.push_to(data_matrix2);
 * 
 *     for(int i = 0; i < params.n_rows; i++){
 *         for(int j = 0; j < params.n_cols; j++){
 *             data_matrix2t(j, i) = data_matrix2(i, j);
 *         }
 *     }
 * 
 *     std::string sdata2 = "output/data22.txt";
 *     std::ofstream d2(sdata2.c_str());
 *     data_matrix2.print(d2, 12 / *width; avoids that two numbers appear as one in the output* /);
 *     std::string sdata2t = "output/gnuplot_inputwaves.txt";
 *     std::ofstream d2t(sdata2t.c_str());
 *     data_matrix2t.print(d2t, 12 / *width; avoids that two numbers appear as one in the output* /);
 * 
 * @endcode
 * 
 * write clusters in separate files
 * 

 * 
 * 
 * @code
 *     std::cout  << std::endl << " ---------- DONE ----------" << std::endl;
 * 
 * #ifdef USE_KMEANS
 * 
 *     std::cout  << std::endl << " ---------- KMEANS ----------" << std::endl;
 * 
 *     std::cout << "n_clusters= " << params.n_clusters <<std::endl;
 * 
 *     step6::Kmeans<double, cublas> cluster_search(data,params.mmeans,params.max_iter);
 *     cluster_search.initial_computation(params.method);
 *     cluster_search.iterate(params.method,params.max_iter);
 * 
 * @endcode
 * 
 * write the calculated means and covariances to the files means.txt and sig.txt respectively
 * 
 * @code
 *     std::string smmeans = "output/mmeans.txt";
 *     std::ofstream j(smmeans.c_str());
 *     params.mmeans.print(j, 12 / *width; avoids that two numbers appear as one in the output* /);
 * 
 *     std::cout  << std::endl << " ---------- DONE ----------" << std::endl;
 * 
 * #endif
 * 
 *     dealii::ParameterHandler prm_debug_handler;
 * 
 *     FullMatrixAccessor means(params.mmeans, true);
 * 
 *     std::cout  << std::endl << " ---------- SHOHAM ----------" << std::endl;
 * 
 *     step7::ShohamEM<T, cublas> Test(data, means); // data, means
 * 
 *     Test.run_Shoham();
 * 
 *     std::cout  << std::endl << " ---------- DONE WITH SHOHAM----------" << std::endl;
 * 
 * #ifdef USE_KMEANS
 * @endcode
 * 
 * write the calculated means and covariances to the files means.txt and sig.txt respectively
 * 
 * @code
 *     std::string smeans = "output/means.txt";
 *     std::ofstream m(smeans.c_str());
 *     params.mmeans.print(m, 12 / *width; avoids that two numbers appear as one in the output* /);
 * #endif
 * 
 * @endcode
 * 
 * write the calculated original waveforms
 * 
 * @code
 *     std::string oforms = "output/oforms.txt";
 *     std::ofstream of(oforms.c_str());
 *     params.check_waves.print(of, 12 / *width; avoids that two numbers appear as one in the output* /);
 * 
 * #ifdef USE_KMEANS
 * @endcode
 * 
 * Write the waveforms to File
 * 
 * @code
 *     dealii::MatrixOut means_out;
 *     std::string means_forms = params.output_folder+"/means-waveforms.txt";
 *     std::ofstream mf(means_forms.c_str());
 *     means_out.build_patches(params.mmeans, "waveforms");
 *     means_out.write_gnuplot(mf);
 * #endif
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Funktionrun_simulation"></a> 
 * <h4>Funktion: run_simulation</h4>
 * 
 * 
 * Pass the user inputs from the main function to the simulation
 * 
 * @code
 * void run_simulation(int / *argc* /, char *argv[])
 * {
 * @endcode
 * 
 * the input parameters and reads from a file named like the program with the extension .prm.
 * The read parameters are written to a log-file.
 * 
 * @code
 *     dealii::ParameterHandler prm_handler;
 * 
 * #ifdef USE_SINGLE_PRECISION
 *     typedef float NumberType;
 * #else
 *     typedef double NumberType;
 * #endif
 * 
 *     typedef step7::SimParams<NumberType> SimParams;
 * 
 *     SimParams::declare(prm_handler);
 * 
 *     std::string prm_filename(argv[0]);
 *     prm_filename += ".prm";
 *     prm_handler.read_input (prm_filename);
 * 
 *     prm_filename += ".log";
 *     std::ofstream log_out_text(prm_filename.c_str());
 *     prm_handler.print_parameters (log_out_text,
 *                                   dealii::ParameterHandler::Text);
 * 
 *     SimParams params;
 * 
 *     params.get(prm_handler);
 *     int DevNo = 0;
 * 
 *     cudaSetDevice(DevNo);           // select CUDA device
 * 
 *     global_data.current_device_id = DevNo;
 * 
 *     global_data.cublanc_gpu_info(); // obtain informations about used GPU
 * 
 * #ifdef USE_SINGLE_PRECISION
 *     std::cout << "Simulation with float:\n"
 *               << "------------------------------------------------" << std::endl;
 * 
 * @endcode
 * 
 * Test with float
 * 
 * @code
 *     step7::MyFancySimulation<float> simu_float(prm_filename);
 * 
 *     simu_float.run_pca(prm_handler);
 * 
 *     simu_float.run_cluster_analysis();
 * 
 *     std::cout << "\nSimulation with float DONE \n"
 *               << "------------------------------------------------\n" << std::endl;
 * #else
 *     std::cout << "Simulation with double:\n"
 *               << "------------------------------------------------" << std::endl;
 * 
 *     fflush(stdout);
 * 
 * @endcode
 * 
 * Test with double
 * 
 * @code
 *     step7::MyFancySimulation<double> simu_double(prm_filename);
 * 
 *     simu_double.run_pca(prm_handler);
 * 
 *     simu_double.run_cluster_analysis();
 * 
 *     std::cout << "\nSimulation with double DONE \n"
 *               << "------------------------------------------------\n" << std::endl;
 * #endif
 * 
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Funktionrun_kerneltest"></a> 
 * <h4>Funktion: run_kerneltest</h4>
 * 
 * 
 * Pass the user inputs from the main function to the simulation
 * 
 * @code
 * void run_kerneltest(int / *argc* /, char *argv[])
 * {
 * 
 *     int DevNo = 0;
 * 
 *     cudaSetDevice(DevNo);
 * #ifdef USE_SINGLE_PRECISION
 *     typedef float NumberType;
 * #else
 *     typedef double NumberType;
 * #endif
 * 
 *     static const int n_data_points = 100000;
 *     static const int n_sampling_points = 64;
 *     static const int n_clusters = 10;
 * 
 *     step7::KernelTest<NumberType, cublas> Test(n_data_points,
 *                                                n_sampling_points,
 *                                                n_clusters);
 * 
 *     Test.run();
 * 
 * }
 * 
 * namespace step7 {
 * 
 * }// namespace step7 END
 * 
 * @endcode
 * 
 * 
 * <a name="Funktionmain"></a> 
 * <h4>Funktion: main</h4>
 * 
 * 
 * start of the simulation
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *     cudaGetDeviceCount(&global_data.n_CUDA_devices);
 *     std::cout
 *             << "N available CUDA devices: "
 *             << global_data.n_CUDA_devices << std::endl;
 * 
 *     run_simulation(argc, argv);
 * 
 * }
 * 
 * 
 * 
 * @endcode
<a name="Results"></a><h1>Results</h1>

The actual not sorted input is depicted in the following plot.
<p>
\image html waveforms3d.png "Input waveforms generated by the generate_waveforms()-function (3 different shapes)" width=15cm
</p>
Now the algorithm should work on this input data.<br /><br />
<b>Initialization phase:</b>
<p>
\htmlimage{progress_0.png, 700, Visualization of the algorithm at the initialization phase}
</p>
<ul>
<li>Upper part:
<ul>
<li>2D plot of the generated waveforms</li>
<li>x-axis: time</li>
<li>y-axis: potential</li>
<li>color indicates the (actual) membership of a certain waveform to a cluster</li>
</ul>
</li>
<li>Middle part:
<ul>
<li>2D plot of the given waveforms transformed into principal components</li>
<li>x-axis: PC1</li>
<li>y-axis: PC2</li>
<li>color indicates the (actual) membership of a certain waveform to a cluster</li>
<li>ellipses represent the (actual) parameters of the clusters (see figure below)</li>
</ul>
</li>
<li>Lower part:
<ul>
<li>bar plot of the $\pi$ values</li>
<li>x-axis: j (number of respective cluster)</li>
<li>y-axis: $\pi$ value</li>
<li>color indicates the (actual) cluster j</li>
</ul></li>
</ul>
<p>
\htmlimage{ellipse_details.png, 700, Cluster parameters indicated by an ellipse}
</p>
<ul>
<li>points are the transformed waveforms $x_i$</li>
<li>blue stars are the means of a ellipse j $\mu_j$</li>
<li>the ellipse j is spanned by its semiaxis (computed of eigenvalues and eigenvectors of $\Sigma_j$) $\Sigma_{j_x}, \Sigma_{j_y}$</li>
<li>ellipse j represents $\Sigma_j$</li>
</ul>
<b>Iteration:</b>
<p>
\htmlimage{progress_1.png, 700, Visualization of the algorithm in iteration phase (optimizing probability with given clusters)}
</p>
<p>
\htmlimage{progress_2.png, 700, Visualization of the algorithm in iteration phase (purging clusters with less points)}
</p>
<p>
\htmlimage{progress_3.png, 700, Visualization of the algorithm in iteration phase (clusters have been reduced to a small subset of high probability clusters)}
</p>
<p>
\htmlimage{progress_4.png, 700, Visualization of the algorithm in iteration phase (final result showing final classifications and parameters)}
</p>
<b>Video:</b>
<p>
\image html video.gif "Visualization of the algorithm" width=15cm
</p>
In order to produce the video after the algorithm was applied, the folloing steps are required:
<ul>
<li>while the algorithm is running several snapshots are generted by the dump() function in the "/plot" directory
<ul>
<li>plots in postscript are generated using gnuplot</li>
<li>postscript files are converted into pdf (postscript removed)</li>
<li>pdf files are rasterized to jpeg files which are copied several times (makes video slower)</li>
</ul>
</li>
<li>to compile all the jpeg files into a video, use the following command:
<pre>ffmpeg -i out\%d.jpeg video.mpg</pre></li>
<li>the dump() function can be called on a abitrary step of the algorithm to show snapshots of interest (or disabled to ensure effiency)</li>
</ul>
<b>Simulation parameters:</b>

The visualized simulation above was executed with the following parameters
<pre>

# Listing of Parameters
# ---------------------
subsection Algorithmic parameter.
  # Maximal number of iterations.
  set Max iterations           = 50

  # cpu or gpu
  set Method                   = cpu

  # Number of openMP threads.
  set Number of openMP threads = 16
end


subsection Dimensions of test problems.
  # Binary logarithm of minimal number of columns of upper triangular matrix
  # R. Allowed range : [3,13]
  set log2(max n cols)   = 3

  # Binary logarithm of minimal number of columns of upper triangular matrix
  # R. Allowed range : [2,12]
  set log2(min n cols)   = 2

  # Repeat the test for a given matrix size this many times.
  set n repetitions      = 1

  # In each test instance the number of rows of R is increased by this factor.
  # Allowed range : [1.1,10]
  set n rows growth rate = 1.2
end


subsection Global parameters
  # double
  set Run Simulation double = false

  # Single precision
  set Run Simulation float  = true
end


subsection NIPALS parameters
  # Maximal number of iterations. If -1 is given termination of NIPALS
  # iterations solely depends on the given tolerance
  set Maximal number of iterations = 1000

  # Number of components
  set Number of components         = 2

  # Once the relative change of an eigenvalue in the inner power iteration is
  # les than the given value iteration stops and the next principal component
  # is computed.
  set Tolerance for eigenvalues    = 1e-5
end


subsection Problem parameter.
  # Cluster centers.
  set Cluster centers                   = (0,0),(3,0),(0,4),(5,5)

  # Initial Cluster centers.
  set Initial centers                   = (0,0),(1,1),(2,2),(3,3)

  # Number of clusters.
  set Number of clusters                = 20

  # Number of data points per cluster.
  set Number of data points per cluster = 1000

  # standard deviation for the normal distribution Allowed range : [0.,10.]
  set Sigma                             = 0.3

  # Corresponds to the number of the principal components one is looking for.
  set dimension                         = 2
end


subsection Program parameter.
  # input file for calculation
  set Input file    = bla.data

  # output folder location
  set Output folder = output
end


subsection Visualization flags.
  # Multiplying Q and R as obtained from the factorization should give
  # something that closely resembles the original A
  set print QR                     = false

  # Computing ||Q_orig - Q||_2 must be done column-wise as the individual
  # columns of Q are reproduced only up to a sign (the algorithm cannot
  # distinguish the direction of a column vector).
  set print Q column errors        = false
  set print Q from factorization   = false

  # Multiplying the Q obtained from the factorization with its own transpose
  # should give a unit matrix.
  set print QtQ from factorization = false
  set print R from factorization   = false

  # The original A is given multiplying the original Q with the original R.
  set print original A             = false
  set print original Q             = false
  set print original R             = false
end


subsection Waveforms parameters
  # Center 1 of Wave 3
  set Center 1 of Wave 3 = 0.6

  # Center 2 of Wave 3
  set Center 2 of Wave 3 = 0.65

  # Center of Wave 1
  set Center of Wave 1   = 0.25

  # Center of Wave 2
  set Center of Wave 2   = 0.85

  # Noise intensity
  set Noise intensity    = 0.2

  # Number of columns
  set Number of columns  = 128

  # Number of Waves
  set Number of waves    = 300,300,300

  # output folder location
  set Output folder      = output

  # Sigma 1 of Wave 3 (Lorenz)
  set Sigma 1 of Wave 3  = 0.05

  # Sigma 2 of Wave 3
  set Sigma 2 of Wave 3  = 0.06

  # Sigma of Wave 1 (standard deviation of gaussian)
  set Sigma of Wave 1    = 0.09
end

</pre>

<b>Memory throughput:</b>

The memory bandwidth and runtimes measurements of the implemented kernels are shown in the next figure. They are obtained on a Tesla C2070 graphic card. Denote that the maximal bandwidth of the hardware is 120 GB/s with 515 GFlop/s (double precision).

<p>
\htmlimage{bandwidth.png, 700, red: processing power of kernels using normal operations\, blue: processing power of kernels using &quot;special function kernels&quot; (SFUs)\, green: memory throughput of kernels.}
</p>

Based on this results it can be obtained that the implemented algorithm is profitable due to the effective use of the graphic hardware bandwidth which is 5-10 times larger than the one of a CPU. The following table shows the runtimes of each kernel.

<table>
<tr><th>Kernel</th><th>Time in ms</th></tr>
<tr><td>Subtract Array</td><td>1.24</td></tr>
<tr><td>Matrix-Matrix Skalar</td><td>10</td></tr>
<tr><td>Determinate</td><td>0.05</td></tr>
<tr><td>Assemble P_ij</td><td>0.96</td></tr>
<tr><td>Assemble u_ij</td><td>0.19</td></tr>
<tr><td>Assemble z_ij</td><td>0.32</td></tr>
<tr><td>Calculate y</td><td>0.38</td></tr>
<tr><td>Sum Array</td><td>0.07</td></tr>
<tr><td>Sum Array2D</td><td>0.41</td></tr>
<tr><td>Matrix-Matrix Multiplikation</td><td>1.3</td></tr>
<tr><td>SMX</td><td>1.26</td></tr>
</table>
 * <a name="PlainProg"></a>
 * <h1> The plain program</h1>
 * 
 * (If you are looking at a locally installed CUDA HPC Praktikum version, then the
 * program can be found at <i>
 *  .. /.. /testsite / /step-7 /step-cu.cc
 * </i>. Otherwise, this is only
 * the path on some remote server.)
 @code

 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef WAVEFORMS_GENERATOR_H
 * #define WAVEFORMS_GENERATOR_H
 * 
 * #include <vector>
 * #include <time.h>
 * #include <iostream>
 * #include <vector>
 * #include <fstream>
 * #include <QDir>
 * #include <QString>
 * #include <QVector>
 * 
 * #include <boost/random/mersenne_twister.hpp>
 * #include <boost/random/normal_distribution.hpp>
 * #include <boost/random/variate_generator.hpp>
 * #include <boost/random/uniform_real.hpp>
 * #include <boost/math/distributions/students_t.hpp>
 * 
@endcode
 <a name="plain-ClassWaveForms"></a>
@code
 * class WaveForms{
 * 
 * public:
 *     WaveForms();
 * 
 * 
@endcode
 <a name="plain-Functiongenerate_waveforms"></a>
@code
 *     template <typename D>
 *     static void generate_waveforms(unsigned int row,
 *                                    unsigned int col,
 *                                    dealii::FullMatrix<D> &waveforms,
 *                                    const std::vector<int> &n_waves,
 *                                    D S,
 *                                    const std::vector<double> &waves_parameter)
 *     {
 *         boost::mt19937 rng;
 *         boost::uniform_real<D> u(-1., 1.);
 *         boost::variate_generator<boost::mt19937&, boost::uniform_real<D> > gen(rng, u);
 * 
 *         boost::mt19937 rng2;
 *         boost::uniform_real<D> v(0., 1.);
 *         boost::variate_generator<boost::mt19937&, boost::uniform_real<D> > gen2(rng2, v);
 * 
 * 
 *         double form;
 *         int form_id;
 *         double X, X1, X2, sig, sig1, sig2, nu, t, Y, noise, l, t0;
 *         double P0,P1,P2,P3,P4,P5,P6;
 * 
 *         P0 = waves_parameter[0];
 *         P1 = waves_parameter[1];
 *         P2 = waves_parameter[2];
 *         P3 = waves_parameter[3];
 *         P4 = waves_parameter[4];
 *         P5 = waves_parameter[5];
 *         P6 = waves_parameter[6];
 * 
 *         for(int j = 0; j < row ; j++){
 * 
 *             form = gen2();
 * 
 *             if(form<(n_waves[0]/(double)row))
 *                 form_id=1;
 *             else{
 *                 if(form > ((n_waves[0]+n_waves[1])/(double)row))
 *                     form_id = 3;
 *                 else
 *                     form_id= 2;
 *             }
 * 
 *             switch(form_id){
 *             case 1:{
 *                 noise = 1+(gen()*S);
 * 
 *                 X = P0*(noise);
 *                 sig = P1*(noise);
 *                 t = 0.1*(noise);
 *                 l = 0.;
 * 
 *                 for(int k = 0; k<col; k++){
 *                     Y = gen()*S;
 *                     l = k/double(col);
 *                     waveforms(j,k) = (1+Y)+(1/t)*((l-X)*exp(-(std::pow(((l-X)/(sig*sqrt(2))),2))));
 *                 }
 *                 break;
 *             }
 *             case 2:{
 * 
 *                 noise = (gen() * S);
 *                 t0 = P2*double(col);
 * 
 *                 for(int k = 0; k<col; k++){
 * 
 *                     t = 12;
 *                     if(k < t0-t/2 || k > t0+t/2)
 *                         waveforms(j,k) = 0.0;
 *                     else
 *                         waveforms(j,k) = (1 + noise) * std::pow(cos(((k - t0 + t * noise) / t) * M_PI), 2);
 *                     waveforms(j,k)   += gen() * S;
 *                 }
 * 
 * 
 *                 break;
 *             }
 *             case 3:{
 * 
 *                 noise = 1+(gen()*S)*0.5;
 * 
 *                 X1 = P3*noise;
 *                 sig1 = P5*noise;
 *                 t = noise;
 *                 X2 = P4*noise;
 *                 sig2 = P6*noise;
 *                 l = 0.;
 * 
 *                 for(int k = 0; k<col; k++){
 *                     Y = gen()*S;
 *                     l = k/double(col);
 *                     waveforms(j,k) = (1+Y)+(1/t)*(((1/(1+std::pow(((l-X1)/sig1),2))) - (1/(1+std::pow(((l-X2)/sig2),2)))));
 *                 }
 *                 break;
 *             }
 *             default:{
 *                 std::cout<<"illegal waveform choosen"<<std::endl;
 *             }
 * 
 *             }
 *         }
 *     }
 * 
@endcode
 <a name="plain-Functionnormalize_waveforms"></a>
@code
 *     template <typename D>
 *     static void normalize_waveforms(unsigned int n_rows, unsigned int n_cols,
 *                                     dealii::FullMatrix<D> &waveforms)
 *     {
 *         double row_sum, row_mean;
 *         for(int i = 0; i < n_rows; i++){
 *             row_sum = 0., row_mean = 1.;
 *             for(int j = 0; j < n_cols; j++){
 *                 row_sum += waveforms(i,j);
 *             }
 *             row_mean = row_sum/(double)n_cols;
 *             for(int j = 0; j < n_cols; j++){
 *                 waveforms(i,j) -= row_mean;
 *             }
 *         }
 *     }
 * 
@endcode
 <a name="plain-Functionoriginal_waveforms"></a>
@code
 *     template <typename D>
 *     static void original_waveforms(unsigned int col, int n_forms,
 *                                    dealii::FullMatrix<D> &check_waves,
 *                                    const std::vector<double> &waves_parameter,
 *                                    const std::vector<int> &n_waves)
 *     {
 *         double X, X1, X2, sig, sig1, sig2, nu, t, noise, l;
 *         double P0,P1,P2,P3,P4,P5,P6;
 * 
 *         P0 = waves_parameter[0];
 *         P1 = waves_parameter[1];
 *         P2 = waves_parameter[2];
 *         P3 = waves_parameter[3];
 *         P4 = waves_parameter[4];
 *         P5 = waves_parameter[5];
 *         P6 = waves_parameter[6];
 * 
 *         int w1, w2, w3;
 * 
 *         w1 = n_waves[0];
 *         w2 = n_waves[1];
 *         w3 = n_waves[2];
 * 
 *         for(int i=0; i<n_forms;i++){
 *             if(w1!= 0 && i<n_forms){
 *                 X = P0;
 *                 sig = P1;
 *                 t = 0.1;
 *                 l = 0.;
 * 
 *                 for(int k = 0; k<col; k++){
 *                     l = k/double(col);
 *                     check_waves(i,k) = (1/t)*((l-X)*exp(-(std::pow(((l-X)/(sig*sqrt(2))),2))));
 *                 }
 *                 w1 = 0;
 *                 i ++;
 *             }
 *             if(w2!=0 && i<n_forms){
 * 
 *                 nu = P2;
 *                 X = (1/(nu*2));
 *                 l= 0.;
 * 
 *                 for(int k = 0; k<col; k++){
 *                     l = k/double(col) + noise;
 *                     if(l > 2*X|| l < 0 ){
 *                         check_waves(i,k) = 0.;
 *                     }else{
 *                         check_waves(i,k) = std::pow(sin(nu*M_PI*l),2);
 *                     }
 *                 }
 *                 w2 =0;
 *                 i++;
 *             }
 *             if(w3!=0 && i<n_forms){
 * 
 *                 X1 = P3;
 *                 sig1 = P5;
 *                 X2 = P4;
 *                 sig2 = P6;
 *                 l = 0.;
 * 
 *                 for(int k = 0; k<col; k++){
 *                     l = k/double(col);
 *                     check_waves(i,k) = (((1/(1+std::pow(((l-X1)/sig1),2))) - (1/(1+std::pow(((l-X2)/sig2),2)))));
 *                 }
 *                 w3 = 0;
 *                 i++;
 *             }
 *         }
 *     }
 * };
 * 
 * 
 * #endif // WAVEFORMS_GENERATOR_H
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef CUDA_KERNEL_STEP_7_CU_H
 * #define CUDA_KERNEL_STEP_7_CU_H
 * 
@endcode
 <a name="plain-Declarationsofthekernelwrapperfunctions"></a>
@code
 * template <typename T>
 * struct Kernel {
 * 
 *     void subtract_array_from_matrix(int n, int p, int j, int g,
 *                                     const T * X,
 *                                     const T * mu_j,
 *                                     T * M);
 * 
 *     void MM_scalar(int n, int p,
 *                    const T * A,
 *                    const T * B,
 *                    T * result);
 * 
 *     T Det(int p, int MS, const T * A);
 * 
 * 
 * 
 *     void assemble_P_ij(T gamma,
 *                        T v,
 *                        int p,int n, int g,
 *                        const T * Delta,
 *                        const T * delta,
 *                        T * P_ij);
 * 
 * 
 *     void assemble_u_ij(T v,
 *                        int p,int n, int g,
 *                        const T * delta,
 *                        T * u_ij);
 * 
 *     void assemble_z_ij(int n, int g,
 *                        const T * P_ij,
 *                        const T * c,
 *                        const T * pi,
 *                        T * z_ij);
 * 
 *     void reduce(unsigned int n,
 *                 const T * c,
 *                 T * result,
 *                 int threads,
 *                 int blocks);
 * 
 *     double calculate_y(int n,
 *                        int g,
 *                        T v,
 *                        T digamma,
 *                        const T * z_ij,
 *                        const T * delta_ij,
 *                        const T * u_ij);
 * 
 *     double sum_array(int n,
 *                      const T * input,
 *                      bool log);
 * 
 *     void sum_matrix_cols(int n, int g,
 *                          const T * input,
 *                          T * output,
 *                          bool log);
 * 
 *     double sum_array(int n,
 *                      int g,
 *                      const T * input,
 *                      bool log);
 * 
 *     void sum_array2D(int n,
 *                      int g,
 *                      const T * z_ij,
 *                      const T * u_ij,
 *                      T * output);
 * 
 *     void MxM(int n,int g,
 *              const T * z_ij,
 *              const T * u_ij,
 *              T * q_ij);
 * 
 *     void SMX(int n, int p,
 *              const T * q_ij,
 *              const T C_j,
 *              const int j,
 *              const T * X,
 *              const T * mu_j,
 *              T * M);
 * 
 * };
 * 
 * #endif // CUDA_KERNEL_STEP_7_CU_H
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #include <step-7/cuda_kernel_wrapper_step-7.cu.h>
 * 
 * #include <math.h>
 * #include <stdio.h>
 * #include <iostream>
 * #include <fstream>
 * 
 * float Giga_bytes = std::pow(1024.,3);
 * 
 * cudaEvent_t beginEvent;
 * cudaEvent_t endEvent;
 * 
 * std::ofstream dt("output/Memory_Bandwidth.txt");
 * 
@endcode
 <a name="plain-Kernel__subtract_array_from_matrix"></a>
@code
 * template <typename T>
 * __global__ void __subtract_array_from_matrix(int n, int p, int j, int g,
 *                                              const T * X,
 *                                              const T * mu_j,
 *                                              T * M)
 * {
 *     int r = blockIdx.y * blockDim.x + threadIdx.x;
 *     int c = blockIdx.x;
 * 
 *     if(r < n && c < p){
 *         M[c*n+r] = X[c*n+r]-mu_j[j+c*g];
 *     }
 * }
 * 
@endcode
 <a name="plain-Kernel__MM_scalar"></a>
@code
 * template <typename T>
 * __global__ void __MM_scalar(int n, int p,
 *                             const T * A,
 *                             const T * B,
 *                             T * result)
 * {
 *     int r = blockIdx.x * blockDim.x + threadIdx.x;
 * 
 *     T sum = 0.;
 * 
 *     if(r < n){
 *         for(int j= 0; j<p; j++){
 *             sum += A[j*n+r]*B[r*p+j];
 *         }
 *         result[r] = sum;
 *     }
 * }
 * 
@endcode
 <a name="plain-Kernel__Det"></a>
@code
 * template <typename T>
 * __global__ void __Det(int p, int MS,
 *                       const T * A,
 *                       T * result)
 * {
 *     T det = 1.;
 * 
 *     for(int i = 0; i<p; i++){
 *         det *= A[i*MS+i];
 *     }
 * 
 *     *result = det;
 * 
 * }
 * 
@endcode
 <a name="plain-Kernel__P_ij"></a>
@code
 * template <typename T>
 * __global__ void __P_ij(T gamma,
 *                        T v,
 *                        int p, int n, int g,
 *                        const T * Delta,
 *                        const T * delta,
 *                        T * P_ij)
 * {
 *     int i = blockIdx.y * blockDim.x + threadIdx.x;
 *     int j = blockIdx.x;
 * 
 * 
 *     if(i < n && j < g){
 *         P_ij[j*n+i] = gamma/Delta[j]/(pow((1.+delta[j*n+i]/v),(v+p)/2.));
 *     }
 * 
 * 
 * }
 * 
@endcode
 <a name="plain-Kernel__u_ij"></a>
@code
 * template <typename T>
 * __global__ void __u_ij(T v,
 *                        int p, int n, int g,
 *                        const T * delta,
 *                        T * u_ij)
 * {
 *     int i = blockIdx.y * blockDim.x + threadIdx.x;
 *     int j = blockIdx.x;
 * 
 * 
 *     if(i < n && j < g){
 *         u_ij[j*n+i] = (p+2)/(delta[j*n+i]+v);
 *     }
 * 
 * }
 * 
@endcode
 <a name="plain-Kernel__z_ij"></a>
@code
 * template <typename T>
 * __global__ void __z_ij(int n, int g,
 *                        const T * P_ij,
 *                        const T * c,
 *                        const T * pi,
 *                        T * z_ij)
 * {
 *     __shared__ T cs[512];
 *     __shared__ T pi_s[32];
 * 
 *     int i = threadIdx.x + blockDim.x * blockIdx.x;
 *     int j = threadIdx.y;
 * 
 * 
 *     if (j == 0){
 *         cs[threadIdx.x] = c[i];
 *     }
 * 
 *     if(threadIdx.x < g && j == 0){
 *         pi_s[threadIdx.x] = pi[threadIdx.x];
 *     }
 * 
 *     __syncthreads();
 * 
 *     if(i < n && j < g){
 *         z_ij[j*n+i] = (P_ij[j*n+i]*pi_s[j])/cs[threadIdx.x];
 *     }
 * 
 * }
 * 
@endcode
 <a name="plain-Functionreduction"></a>
@code
 * template<typename T>
 * __device__ void reduction(volatile T * sum, int tid)
 * {
 *     for(int k = blockDim.x/2; k>0; k/=2){
 *         if(tid < k){
 *             sum[tid] += sum[tid + k];
 *         }
 *         __syncthreads();
 *     }
 * }
 * 
@endcode
 <a name="plain-Kernel__sum_array"></a>
@code
 * template <typename T>
 * __global__ void __sum_array(int n,
 *                             const T * input,
 *                             T * result)
 * {
 *     __shared__ T sum[512];
 * 
 *     sum[threadIdx.x] = 0.;
 *     __syncthreads();
 * 
 * 
 *     int i = threadIdx.x + blockDim.x * blockIdx.x;
 * 
 *     if(i>=n)
 *         return;
 * 
 *     sum[threadIdx.x] += input[i];
 *     __syncthreads();
 * 
 *     reduction(sum, threadIdx.x);
 * 
 *     if(threadIdx.x == 0){
 *         result[blockIdx.x]= sum[0];
 *     }
 * 
 * }
 * 
@endcode
 <a name="plain-Kernel__sum_log_array"></a>
@code
 * template <typename T>
 * __global__ void __sum_log_array(int n,
 *                                 const T * input,
 *                                 T * result)
 * {
 *     __shared__ T sum[512];
 * 
 *     sum[threadIdx.x] = 0.;
 *     __syncthreads();
 * 
 * 
 *     int i = threadIdx.x + blockDim.x * blockIdx.x;
 * 
 *     if(i>=n)
 *         return;
 * 
 *     sum[threadIdx.x] += log(input[i]);
 *     __syncthreads();
 * 
 *     reduction(sum, threadIdx.x);
 * 
 *     if(threadIdx.x == 0){
 *         result[blockIdx.x]= sum[0];
 *     }
 * 
 * }
 * 
@endcode
 <a name="plain-Kernel__sum_array2D"></a>
@code
 * template <typename T>
 * __global__ void __sum_array2D(int n, int g,
 *                               const T * z_ij,
 *                               const T * u_ij,
 *                               T * output)
 * {
 *     __shared__ T sum[512];
 * 
 *     sum[threadIdx.x] = 0.;
 *     __syncthreads();
 * 
 * 
 *     int i = threadIdx.x + blockDim.x * blockIdx.x;
 *     int j = blockIdx.y;
 * 
 *     if(i>=n)
 *         return;
 * 
 *     int global_index = j*n+i;
 *     sum[threadIdx.x] += z_ij[global_index]* u_ij[global_index];
 *     __syncthreads();
 * 
 *     reduction(sum, threadIdx.x);
 * 
 *     if(threadIdx.x == 0){
 *         output[blockIdx.y*gridDim.x+blockIdx.x]= sum[0];
 *     }
 * 
 * }
 * 
@endcode
 <a name="plain-Kernel__y"></a>
@code
 * template <typename T>
 * __global__ void __y(int n, int g, T v, T digamma,
 *                     const T * z_ij,
 *                     const T * delta_ij,
 *                     const T * u_ij,
 *                     T * results)
 * {
 *     __shared__ T sum[256];
 * 
 *     sum[threadIdx.x] = 0.;
 *     __syncthreads();
 * 
 * 
 *     int i = threadIdx.x + blockDim.x * blockIdx.x;
 * 
 *     if(i>=n)
 *         return;
 * 
 *     for(int j = 0; j<g; j++){
 *         int global_index = j*n+i;
 *         sum[threadIdx.x] += z_ij[global_index]*(digamma+log(2/(delta_ij[global_index]+v))-u_ij[global_index]);
 *     }
 *     __syncthreads();
 * 
 *     reduction(sum, threadIdx.x);
 * 
 *     if(threadIdx.x == 0){
 *         results[blockIdx.x]= sum[0];
 *     }
 * 
 * }
 * 
@endcode
 <a name="plain-Kernel__MxM"></a>
@code
 * template <typename T>
 * __global__ void __MxM(int n, int g,
 *                       const T * z_ij,
 *                       const T * u_ij,
 *                       T * q_ij)
 * {
 *     int i = blockIdx.y * blockDim.x + threadIdx.x;
 *     int j = blockIdx.x;
 * 
 * 
 *     if(i < n && j < g){
 *         q_ij[j*n+i] = z_ij[j*n+i]*u_ij[j*n+i];
 *     }
 * }
 * 
@endcode
 <a name="plain-Kernel__SMX"></a>
@code
 * template <typename T>
 * __global__ void __SMX(int n, int p,
 *                       const T * q_ij,
 *                       const T inv_C_j,
 *                       const int j,
 *                       const T * X,
 *                       const T * mu_j,
 *                       T * M)
 * {
 *     int r = blockIdx.y * blockDim.x + threadIdx.x;
 *     int c = blockIdx.x;
 * 
 * 
 *     if(r < n && c < p){
 *         M[c*n+r] = (X[c*n+r]-mu_j[c]) * sqrt(q_ij[j*n+r] * inv_C_j);
 *     }
 * }
 * 
@endcode
 <a name="plain-Functionsubtract_array_from_matrix"></a>
@code
 * template<typename T>
 * void Kernel<T>::subtract_array_from_matrix(int n, int p, int j, int g,
 *                                            const T * X,
 *                                            const T * mu_j,
 *                                            T * M)
 * {
 *     int block_size = 512;
 *     int gridDim_y = (n+block_size-1)/block_size;
 * 
 *     dim3 num_blocks(p, gridDim_y);
 * 
 *     double time_sum = 0;
 *     int n_operations = 1;
 *     int n_memory_access = 3;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * 
 *     dt <<"Kernel name"
 *       <<"\t"
 *      <<"Time in Milliseconds"
 *     <<"\t\t"
 *     <<"bandwidth in GB/s"
 *     <<"\t\t"
 *     <<"Gflops/Second"
 *     << std::endl;
 *     dt << "------------------------------------------------------------------------------------------------------------"
 *        <<std::endl;
 * 
 *     for(int s = 0; s < N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         __subtract_array_from_matrix<<<num_blocks, block_size>>>(n, p, j, g,
 *                                                                  X,
 *                                                                  mu_j,
 *                                                                  M);
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 *     }
 *     time_sum /= N_REPS;
 *     dt << "subtract_kernel:::::"
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*p*sizeof(T))*n_memory_access)/(time_sum/1000))/ Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*g*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 * }
 * 
@endcode
 <a name="plain-FunctionMM_scalar"></a>
@code
 * template<typename T>
 * void Kernel<T>::MM_scalar(int n, int p,
 *                           const T * A,
 *                           const T * B,
 *                           T * output)
 * {
 *     int block_size = 512;
 *     int num_blocks = (n+block_size-1)/block_size;
 * 
 *     double time_sum = 0;
 *     long n_operations = 2*p;
 *     int n_memory_access = 2*p+1;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 *         __MM_scalar<<<num_blocks, block_size>>>(n, p,
 *                                                 A,
 *                                                 B,
 *                                                 output);
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 *     }
 * 
 *     time_sum /= N_REPS;
 *     dt << "MM_scalar:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*p*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*p)*n_operations)/(time_sum/1000)/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * }
 * 
@endcode
 <a name="plain-FunctionDet"></a>
@code
 * template<typename T>
 * T
 * Kernel<T>::Det(int p, int MS, const T * A)
 * {
 *     int block_size = 1;
 *     int num_blocks = 1;
 * 
 *     double time_sum = 0;
 *     int n_operations = p;
 *     int n_memory_access = p+1;
 * 
 * 
 *     T det;
 *     T * det_d;
 *     cudaMalloc(&det_d, sizeof(T));
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         __Det<<<num_blocks, block_size>>>(p, MS, A, det_d);
 * 
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 *     }
 *     time_sum /= N_REPS;
 *     dt << "Det_Kernel:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((p*p*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((p*p*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 *     cudaThreadSynchronize();
 * 
 *     cudaMemcpy(&det, det_d,
 *                sizeof(T), cudaMemcpyDeviceToHost);
 * 
 *     cudaFree(det_d);
 * 
 *     return det;
 * }
 * 
@endcode
 <a name="plain-Functionassemble_P_ij"></a>
@code
 * template<typename T>
 * void
 * Kernel<T>::assemble_P_ij(T gamma,
 *                          T v,
 *                          int p, int n, int g,
 *                          const T * Delta,
 *                          const T * delta,
 *                          T * P_ij)
 * {
 *     int block_size = 512;
 *     int gridDim_y = (n+block_size-1)/block_size;
 * 
 *     dim3 num_blocks(g, gridDim_y);
 * 
 * 
 *     double time_sum = 0;
 *     int n_operations = 4;
 *     int n_memory_access = 3;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 * 
 *         __P_ij<<<num_blocks, block_size>>>(gamma,
 *                                            v,
 *                                            p, n, g,
 *                                            Delta,
 *                                            delta,
 *                                            P_ij);
 * 
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * 
 *     }
 *     time_sum /= N_REPS;
 *     dt << "assemble_P_ij:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*g*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*g*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 * }
 * 
@endcode
 <a name="plain-Functionassemble_u_ij"></a>
@code
 * template<typename T>
 * void
 * Kernel<T>::assemble_u_ij(T v,
 *                          int p, int n, int g,
 *                          const T * delta,
 *                          T * u_ij)
 * {
 *     int block_size = 512;
 *     int gridDim_y = (n+block_size-1)/block_size;
 * 
 *     dim3 num_blocks(g, gridDim_y);
 * 
 * 
 *     double time_sum = 0;
 *     int n_operations = 2;
 *     int n_memory_access = 2;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * 
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         __u_ij<<<num_blocks, block_size>>>(v,
 *                                            p, n, g,
 *                                            delta,
 *                                            u_ij);
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * 
 *     }
 *     time_sum /= N_REPS;
 *     dt << "assemble_u_ij:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*g*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*g*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 * }
 * 
@endcode
 <a name="plain-Functionassemble_z_ij"></a>
@code
 * template<typename T>
 * void
 * Kernel<T>::assemble_z_ij(int n, int g,
 *                          const T * P_ij,
 *                          const T * c,
 *                          const T * pi,
 *                          T * z_ij)
 * {
 *     int blockDim_x = max(16/g, 1) * 32;
 *     dim3 block_size(blockDim_x,g);
 *     int num_blocks = (n+blockDim_x-1)/blockDim_x;
 * 
 * 
 *     double time_sum = 0;
 *     int n_operations = 2;
 *     int n_memory_access = 2;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * 
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 *         __z_ij<<<num_blocks, block_size>>>(n, g,
 *                                            P_ij,
 *                                            c,
 *                                            pi,
 *                                            z_ij);
 * 
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * 
 *     }
 *     time_sum /= N_REPS;
 *     dt << "assemble_z_ij:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << ((((n*g*sizeof(T))*n_memory_access)+n+g)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*g*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 * }
 * 
@endcode
 <a name="plain-Functioncalculate_y"></a>
@code
 * template<typename T>
 * double
 * Kernel<T>::calculate_y(int n, int g, T v, T digamma,
 *                        const T * z_ij,
 *                        const T * delta_ij,
 *                        const T * u_ij)
 * {
 *     int block_size = 256;
 *     int num_blocks = (n+block_size-1)/block_size;
 * 
 *     T y;
 * 
 *     T * results_h = new T[num_blocks];
 *     T * results_d;
 *     cudaMalloc(&results_d, num_blocks*sizeof(T));
 * 
 *     double time_sum = 0;
 *     int n_operations = 6*g + sqrt(block_size/2)+1;
 *     int n_memory_access = 4;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * 
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         __y<<<num_blocks, block_size>>>(n, g, v, digamma,
 *                                         z_ij,
 *                                         delta_ij,
 *                                         u_ij,
 *                                         results_d);
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * 
 *     }
 *     time_sum /= N_REPS;
 *     dt << "calculate_y:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*g*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*g*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 *     cudaThreadSynchronize();
 * 
 *     cudaMemcpy(results_h, results_d,
 *                num_blocks*sizeof(T), cudaMemcpyDeviceToHost);
 * 
 *     for(int i=0; i<num_blocks; i++){
 *         y += results_h[i];
 *     }
 * 
 *     return y;
 * }
 * 
@endcode
 <a name="plain-Functionsum_array"></a>
@code
 * template<typename T>
 * double
 * Kernel<T>::sum_array(int n,
 *                      const T * input,
 *                      bool log)
 * {
 *     double result = 0.;
 *     int block_size = 512;
 *     int num_blocks = (n+block_size-1)/block_size;
 * 
 * 
 *     T * results_h = new T[num_blocks];
 *     T * results_d;
 *     cudaMalloc(&results_d, num_blocks*sizeof(T));
 * 
 *     double time_sum = 0;
 *     int n_operations = sqrt(block_size/2)+1+1;
 *     int n_memory_access = 2;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * 
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         if(log){
 *             __sum_log_array<<<num_blocks, block_size>>>(n,input,results_d);
 *         }else{
 *             __sum_array<<<num_blocks, block_size>>>(n,input,results_d);
 *         }
 * 
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * 
 *     }
 *     time_sum /= N_REPS;
 *     dt << "sum_array:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 *     cudaThreadSynchronize();
 * 
 *     cudaMemcpy(results_h, results_d,
 *                num_blocks*sizeof(T), cudaMemcpyDeviceToHost);
 * 
 *     for(int i=0; i<num_blocks; i++){
 *         result += results_h[i];
 *     }
 * 
 *     return result;
 * }
 * 
@endcode
 <a name="plain-Functionsum_matrix_cols"></a>
@code
 * template<typename T>
 * void
 * Kernel<T>::sum_matrix_cols(int n, int g,
 *                            const T * input,
 *                            T * output,
 *                            bool log)
 * {
 *     int block_size = 512;
 *     int num_blocks = (n+block_size-1)/block_size;
 * 
 * 
 *     T * results_h = new T[num_blocks*g];
 *     T * results_d;
 *     cudaMalloc(&results_d, num_blocks*g*sizeof(T));
 * 
 *     for(int j = 0;j < g; j++){
 *         if(log){
 *             __sum_log_array<<<num_blocks, block_size>>>(n,input+(n*j),results_d+(j*num_blocks));
 *         }else{
 *             __sum_array<<<num_blocks, block_size>>>(n,input+(n*j),results_d+(j*num_blocks));
 *         }
 *     }
 *     cudaThreadSynchronize();
 * 
 *     cudaMemcpy(results_h, results_d,
 *                num_blocks*g*sizeof(T), cudaMemcpyDeviceToHost);
 * 
 *     for(int j = 0; j<g; j++){
 *         for(int i=0; i<num_blocks; i++){
 *             output[j] += results_h[i+j*num_blocks];
 *         }
 *     }
 * 
 * }
 * 
@endcode
 <a name="plain-Functionsum_array2D"></a>
@code
 * template<typename T>
 * void
 * Kernel<T>::sum_array2D(int n, int g,
 *                        const T * z_ij,
 *                        const T * u_ij,
 *                        T * output)
 * {
 *     int block_size = 512;
 *     int gridDim_x = ((n+block_size-1)/block_size);
 *     dim3 num_blocks (gridDim_x, g);
 * 
 * 
 *     T * results_h = new T[gridDim_x*g];
 *     T * results_d;
 *     cudaMalloc(&results_d, gridDim_x*g*sizeof(T));
 * 
 *     double time_sum = 0;
 *     int n_operations = 2 + sqrt(block_size/2)+1;
 *     int n_memory_access = 2;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * 
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         __sum_array2D<<<num_blocks, block_size>>>(n,g,z_ij, u_ij,results_d);
 * 
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * 
 *     }
 *     time_sum /= N_REPS;
 *     dt << "sum_array2D:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*g*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*g*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 *     cudaThreadSynchronize();
 * 
 *     cudaMemcpy(results_h, results_d,
 *                gridDim_x*g*sizeof(T), cudaMemcpyDeviceToHost);
 * 
 *     for(int i=0; i<g; i++){
 *         for(int j= 0; j<gridDim_x; j++){
 *             output[i] += results_h[i*gridDim_x+j];
 *         }
 *     }
 * 
 * }
 * 
@endcode
 <a name="plain-FunctionMxM"></a>
@code
 * template<typename T>
 * void
 * Kernel<T>::MxM(int n, int g,
 *                const T * z_ij,
 *                const T * u_ij,
 *                T * q_ij)
 * {
 *     int block_size = 512;
 *     int gridDim_y = (n+block_size-1)/block_size;
 * 
 *     dim3 num_blocks(g, gridDim_y);
 * 
 * 
 *     double time_sum = 0;
 *     int n_operations = 1;
 *     int n_memory_access = 3;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * 
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         __MxM<<<num_blocks, block_size>>>(n, g,
 *                                           z_ij,
 *                                           u_ij,
 *                                           q_ij);
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * 
 *     }
 *     time_sum /= N_REPS;
 *     dt << "MxM_Kernel:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*g*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*g*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 * }
 * 
@endcode
 <a name="plain-FunctionSMX"></a>
@code
 * template<typename T>
 * void
 * Kernel<T>::SMX(int n, int p,
 *                const T * q_ij,
 *                const T C_j,
 *                const int j,
 *                const T * X,
 *                const T * mu_j,
 *                T * M)
 * {
 *     int block_size = 512;
 *     int gridDim_y = (n+block_size-1)/block_size;
 * 
 *     dim3 num_blocks(p, gridDim_y);
 * 
 * 
 *     double time_sum = 0;
 *     int n_operations = 4;
 *     int n_memory_access = 5;
 * 
 *     cudaEventCreate(&beginEvent);
 *     cudaEventCreate(&endEvent);
 * 
 * 
 *     for(int s = 0; s<N_REPS; s++){
 * 
 *         cudaEventRecord(beginEvent,0);
 * 
 *         __SMX<<<num_blocks, block_size>>>(n, p,
 *                                           q_ij,
 *                                           T(1./C_j),
 *                                           j,
 *                                           X,
 *                                           mu_j,
 *                                           M);
 *         cudaThreadSynchronize();
 * 
 *         cudaEventRecord(endEvent,0);
 * 
 *         cudaEventSynchronize(endEvent);
 * 
 *         float timeValue;
 *         cudaEventElapsedTime(&timeValue, beginEvent, endEvent);
 *         time_sum += timeValue;
 * 
 * 
 *     }
 *     time_sum /= N_REPS;
 *     dt << "SMX_Kernel:::::  "
 *        <<"\t"
 *       << "Time: "
 *       << time_sum
 *       << "\t\t"
 *       << "Memory bandwidth: "
 *       << (((n*p*sizeof(T))*n_memory_access)/(time_sum/1000))/Giga_bytes
 *       << "\t\t"
 *       << "Gflops: "
 *       << ((n*p*n_operations)/(time_sum/1000))/Giga_bytes
 *       << std::endl;
 * 
 *     cudaEventDestroy(beginEvent);
 *     cudaEventDestroy(endEvent);
 * 
 * }
 * 
 * template class Kernel<float>;
 * 
 * #ifndef USE_SINGLE_PRECISION
 * template class Kernel<double>;
 * #endif
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef CUDADriver_STEP_7_H
 * #define CUDADriver_STEP_7_H
 * 
 * #include <deal.II/base/parameter_handler.h>
 * 
 * #include <lac/blas++.h>
 * 
 * #ifdef USE_KMEANS
 * #include <step-6/cuda_driver_step-6.h>
 * #include <step-6/cuda_driver_step-6.hh>
 * #include <step-6/Kmeans.h>
 * #endif
 * 
 * #include <step-7/cuda_kernel_wrapper_step-7.cu.h>
 * 
 * namespace step7 {
 * 
@endcode
 <a name="plain-ClassShohamEM"></a>
@code
 * template <typename T, typename blas>
 * class ShohamEM : public Kernel<T> {
 * public:
 * public:
 * 
 *     typedef typename blas_pp<T, blas>::blas_wrapper_type  BW;
 *     typedef typename blas_pp<T, blas>::FullMatrixAccessor FullMatrixAccessor;
 *     typedef typename blas_pp<T, blas>::Matrix             Matrix;
 *     typedef typename blas_pp<T, blas>::SubMatrix          SubMatrix;
 *     typedef typename blas_pp<T, blas>::MatrixSubCol       MatrixSubCol;
 *     typedef typename blas_pp<T, blas>::MatrixSubRow       MatrixSubRow;
 * 
 *     typedef typename blas_pp<T, blas>::SubColVector SubColVector;
 *     typedef typename blas_pp<T, blas>::SubRowVector SubRowVector;
 * 
 *     typedef typename blas_pp<T, blas>::Vector             Vector;
 *     typedef          transpose<Matrix>                    tr;
 *     typedef          transpose<SubMatrix>                 trSub;
 * 
 *     Matrix mu;
 *     int g;
 *     int g_min;
 *     std::vector<T> pi;
 *     std::vector<Matrix> Sigma;
 *     T v;
 *     double L;
 *     double L_max;
 *     double N;
 *     int p;
 *     Matrix u_ij;
 *     Matrix z_ij;
 *     Matrix P_ij;
 * 
 *     Matrix data;
 *     Matrix mean;
 *     int n_data_points;
 *     Matrix delta_ij;
 *     std::vector<Matrix> Sigma_inv;
 *     T gamma;
 *     Vector Delta;
 *     std::vector<Matrix> Sigma_dec;
 *     Vector pi_d;
 *     Matrix M;
 *     int iter;
 * 
 *     std::vector<T> pi_opt;
 *     Matrix mu_opt;
 *     std::vector<Matrix> Sigma_opt;
 *     T v_opt;
 *     Matrix u_ij_opt;
 *     Matrix z_ij_opt;
 *     Matrix P_ij_opt;
 * 
 *     int it;
 * 
 *     ShohamEM(FullMatrixAccessor &ddata, FullMatrixAccessor &mmeans);
 * 
 *     void dump();
 * 
 *     void run_Shoham();
 *     void initialize_Shoham();
 *     void e_step_Shoham();
 *     void m_step_Shoham();
 *     void purge_step_Shoham();
 * 
 *     ~ShohamEM();
 * 
 * };
 * 
@endcode
 <a name="plain-ClassKernelTest"></a>
@code
 * template <typename T, typename blas>
 * class KernelTest : public Kernel<T> {
 * public:
 *     typedef typename blas_pp<T, blas>::blas_wrapper_type  BW;
 *     typedef typename blas_pp<T, blas>::FullMatrixAccessor FullMatrixAccessor;
 *     typedef typename blas_pp<T, blas>::Matrix             Matrix;
 *     typedef typename blas_pp<T, blas>::SubMatrix          SubMatrix;
 *     typedef typename blas_pp<T, blas>::MatrixSubCol       MatrixSubCol;
 *     typedef typename blas_pp<T, blas>::MatrixSubRow       MatrixSubRow;
 *     typedef typename blas_pp<T, blas>::Vector             Vector;
 *     typedef          transpose<Matrix>              tr;
 *     typedef          transpose<SubMatrix>              trSub;
 * 
 * 
 *     KernelTest(int n_data_points, int dim, int n_clusters);
 * 
 *     int n_data_points, n_clusters, dim;
 * 
 *     void run();
 * 
 *     ~KernelTest();
 * };
 * 
 * } // namespace step7 END
 * 
 * #endif // CUDADriver_STEP_7_H
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef CUDA_DRIVER_STEP_7_HH
 * #define CUDA_DRIVER_STEP_7_HH
 * #include <cuda_driver_step-7.h>
 * 
 * #include <cuda_kernel_wrapper_step-7.cu.h>
 * 
 * #include <cuda_kernel_wrapper_step-1.cu.h>
 * 
 * #include <base/CUDATimer.h>
 * 
 * #include <deal.II/base/timer.h>
 * 
 * #include <cmath>
 * #include <limits>
 * 
 * #include <lac/cublas_wrapper.hh>
 * 
 * #include <boost/math/special_functions/digamma.hpp>
 * #include <boost/math/special_functions/erf.hpp>
 * #include <boost/math/special_functions/beta.hpp>
 * 
 * #include <deal.II/numerics/histogram.h>
 * 
 * #include <iostream>
 * #include <fstream>
 * 
 * #include <QTime>
 * 
 * #include <deal.II/base/convergence_table.h>
 * #include <deal.II/base/table_handler.h>
 * #include <deal.II/base/timer.h>
 * 
 * #include <fstream>
 * 
@endcode
 <a name="plain-StructRightMSolve"></a>
@code
 * template <typename T, typename blas>
 * struct RightMSolve{
 * 
 *     typedef bw_types::SubMatrixView<T, blas> SubMatrix;
 * 
 *     const SubMatrix &l;
 *     SubMatrix &r;
 * 
 *     RightMSolve(const SubMatrix & _l, SubMatrix & _r):l(_l),r(_r){}
 * 
@endcode
 <a name="plain-Operator"></a>
@code
 *     RightMSolve<T, blas> & operator = (const SubMatrix & rhs)
 *     {
 *         r = rhs.matrix();
 * 
 *         std::cout << "? X : AX = B " << std::endl;
 *         char side = 'L', uplo = 'U', transa = 'N', diag = 'N';
 *         int m = rhs.r_end() - rhs.r_begin();
 *         int n = rhs.c_end() - rhs.c_begin();
 * 
 *         T alpha = 1.0;
 * 
 *         const T * const A = l.val();
 *         int lda = l.leading_dim();
 * 
 *         T *  B = r.val();
 *         int ldb = r.leading_dim();
 *         blas::trsm(  side, uplo, transa, diag,
 *                      m, n, alpha,
 *                      A, lda, B, ldb     );
 * 
 *         return *this;
 *     }
 * 
 * };
 * 
@endcode
 <a name="plain-StructRightMSolveTr"></a>
@code
 * template <typename T, typename blas>
 * struct RightMSolveTr{
 * 
 *     typedef bw_types::SubMatrixView<T, blas> SubMatrix;
 *     typedef          transpose<SubMatrix>           tr;
 * 
 *     const SubMatrix &l;
 *     SubMatrix &r;
 * 
 *     RightMSolveTr(const tr & _l, SubMatrix & _r):l(_l.A),r(_r){}
 * 
@endcode
 <a name="plain-Operator"></a>
@code
 *     RightMSolveTr<T, blas> & operator = (const SubMatrix & rhs)
 * 
 *     {
 * 
 *         r = rhs.matrix();
 * 
 *         std::cout << "? X : A^TX = B " << std::endl;
 *         char side = 'L', uplo = 'U', transa = 'T', diag = 'N';
 *         int m = rhs.r_end() - rhs.r_begin();
 *         int n = rhs.c_end() - rhs.c_begin();
 * 
 *         T alpha = 1.0;
 * 
 *         const T * const A = l.val();
 *         int lda = l.leading_dim();
 * 
 *         T * B = r.val();
 *         int ldb = r.leading_dim();
 *         cublas::trsm(  side, uplo, transa, diag,
 *                        m, n, alpha,
 *                        A, lda, B, ldb     );
 * 
 *         return *this;
 *     }
 * 
 * };
 * 
@endcode
 <a name="plain-Operator"></a>
@code
 * template <typename T, typename blas>
 * inline RightMSolve<T, blas>  operator * (const bw_types::SubMatrixView<T, blas> &l, bw_types::SubMatrixView<T, blas> & r){
 *     return RightMSolve<T, blas>(l,r);
 * }
 * 
@endcode
 <a name="plain-Operator"></a>
@code
 * template <typename T, typename blas>
 * inline RightMSolveTr<T, blas>  operator * (const typename RightMSolveTr<T, blas>::tr  &l,
 *                                            bw_types::SubMatrixView<T, blas> & r){
 *     return RightMSolveTr<T, blas>(l,r);
 * }
 * 
@endcode
 <a name="plain-ConstructorShohamEM"></a>
@code
 * template<typename T, typename blas>
 * step7::ShohamEM<T, blas>::ShohamEM(FullMatrixAccessor &ddata, FullMatrixAccessor &mmeans)
 *     :
 *       n_data_points(ddata.n_rows()), p(ddata.n_cols()),
 *       g(mmeans.n_rows()), data(ddata.n_rows(), ddata.n_cols())
 * {
 *     BW::Init();
 *     data = ddata;
 *     mean = mmeans;
 * }
 * 
@endcode
 <a name="plain-Functiondump"></a>
@code
 * template<typename T, typename blas>
 * void step7::ShohamEM<T,blas>::dump()
 * {
 * 
 * 
 *     std::fstream em;
 *     em.open("output/em.dat", std::ios::out);
 * 
 *     for(int i = 0; i < mu.n_rows(); i++){
 * 
 *         for(int j = 0; j < mu.n_cols(); j++)
 *             em << mu(i, j) << " ";
 * 
 *         em << "\n";
 * 
 *     }
 * 
 *     em.close();
 * 
 * 
 *     std::fstream piout;
 *     piout.open("output/pi.dat", std::ios::out);
 * 
 *     for(int i = 0; i < pi.size(); i++)
 *         piout << pi[i] << "\n";
 * 
 *     piout.close();
 * 
 * 
 *     std::fstream ellipse;
 *     ellipse.open("output/ellipses.dat", std::ios::out);
 * 
 *     for(int k = 0; k < g; k++){
 * 
 *         double d;
 *         double eig[2], eigv[2];
 *         int perm[2];
 * 
 *         d = (Sigma[k](0, 0) - Sigma[k](1, 1)) / 2;
 * 
 *         eig[0] = (Sigma[k](0, 0) + Sigma[k](1, 1)) / 2 + sqrt(d * d + Sigma[k](0, 1) * Sigma[k](1, 0));
 *         eig[1] = (Sigma[k](0, 0) + Sigma[k](1, 1)) / 2 - sqrt(d * d + Sigma[k](0, 1) * Sigma[k](1, 0));
 * 
 *         if(eig[0] > eig[1] || true){
 *             perm[0] = 0;
 *             perm[1] = 1;
 *         }
 *         else{
 *             perm[0] = 1;
 *             perm[1] = 0;
 *         }
 * 
 *         std::cout << "\nEig: " << eig[perm[0]] << "\n";
 * 
 *         eigv[0] = eig[perm[0]] - Sigma[k](1, 1);
 *         eigv[1] = Sigma[k](1, 0);
 * 
 *         std::cout << "Eigv: " << eigv[0] / sqrt(eigv[0] * eigv[0] + eigv[1] * eigv[1]) << "..." << perm[0] << "\n";
 * 
 *         ellipse << mu(k, 0) << " " << mu(k, 1) << " " << / *sqrt(* /sqrt(sqrt(eig[perm[0]]))/ *)* / << " " << sqrt(sqrt(sqrt(eig[perm[1]])))
 *                 << " " <<  acos(eigv[0] / sqrt(eigv[0] * eigv[0] + eigv[1] * eigv[1])) * 180 / M_PI << "\n";
 * 
 *     }
 * 
 *     ellipse.close();
 * 
 * 
 *     QString plot_waves, plot_pca, plot_ellipses, plot_hist, plot_process;
 * 
 *     FILE *gp = popen("/usr/local/lib/FEM/gnuplot/bin/gnuplot", "w");//open a pipe to gnuplot
 * 
 *     if(iter == 0){
 *         plot_waves += "!rm plot/ *\n";
 *     }
 * 
 *     plot_waves += "set terminal postscript landscape enhanced color solid linewidth 1.0 'Helvetica' 15\n"
 *             "set output 'plot/out" + QString::number(iter) + ".ps'\n"
 *             "set multiplot\n"
 *             "set xrange [0:128]\n"
 *             "set yrange [-1.5:1.5]\n"
 *             "set size 1,0.3\n"
 *             "set origin 0,0.8\n"
 *             "set key off\n";
 * 
 *     plot_pca += "reset\n"
 *             "set key off\n"
 *             "set size 1,0.7\n"
 *             "set origin 0,0.15\n"
 *             "set style data points\n"
 *             "set yrange []\n"
 *             "set xrange []\n";
 * 
 *     plot_hist = "set origin 0,-0.1\n"
 *             "set size 1,0.3\n"
 *             "set boxwidth 0.8\n"
 *             "set style fill solid\n"
 *             "set yrange [0:0.3]\n"
 *             "set xrange [-0.5:0.5]\n";
 * 
 *     for(int i = 1; i <= n_data_points; i+=1){
 * 
 *         int max = 0;
 * 
 *         for(int j = 1; j < z_ij.n_cols(); j++)
 *             if(z_ij(i-1, j) > z_ij(i-1, max))
 *                 max = j;
 * 
 *         int outlier = 1;
 *         if(u_ij(i-1, max) < 0.75) outlier = 4;
 *         if(u_ij(i-1, max) < 0.4) outlier = 6;
 * 
 *         outlier = 1;
 * 
 *         if(i == 1){
 *             plot_waves += "plot";
 *             plot_pca += "plot";
 *         }
 *         else{
 *             plot_waves += ",";
 *             plot_pca += ",";
 *         }
 * 
 *         plot_waves += " 'output/gnuplot_inputwaves.txt' u " + QString::number(i) +
 *                 "with lines lt 1 lc " + QString::number(max);
 * 
 *         plot_pca += " 'output/pca_out.txt' every ::" + QString::number(i - 1) + "::" + QString::number(i - 1) +
 *                 " lt " + QString::number(outlier) + " lc " + QString::number(max);
 * 
 *     }
 * 
 *     for(int i = 0; i < g; i++){
 * 
 *         plot_ellipses += ", 'output/ellipses.dat' every ::" + QString::number(i) + "::" + QString::number(i) +
 *                 " lc " + QString::number(i) + " with ellipses";
 * 
 *         if(i == 0)
 *             plot_hist += "plot";
 *         else
 *             plot_hist += ",";
 * 
 *         plot_hist += " 'output/pi.dat' every ::" + QString::number(i) + "::" + QString::number(i) +
 *                 " lc " + QString::number(i) + " with histograms";
 * 
 *     }
 * 
 *     plot_waves += "\n";
 *     fprintf(gp, plot_waves.toStdString().c_str());
 *     fflush(gp);
 * 
 *     fprintf(gp, plot_pca.toStdString().c_str());
 *     fflush(gp);
 * 
 *     fprintf(gp, ", 'output/em.dat' lt 3");
 * 
 *     plot_ellipses += "\n";
 *     fprintf(gp, plot_ellipses.toStdString().c_str());
 *     fflush(gp);
 * 
 *     plot_hist += "\n";
 *     fprintf(gp, plot_hist.toStdString().c_str());
 *     fflush(gp);
 * 
 *     plot_process = "!ps2pdf plot/out" + QString::number(iter) + ".ps plot/out" + QString::number(iter) + ".pdf\n"
 *             "!rm plot/out" + QString::number(iter) + ".ps\n"
 *             "!convert -density 300 plot/out" + QString::number(iter) + ".pdf plot/out" + QString::number(iter) + ".jpeg\n";
 * 
 *     int frames_per_image = 30; // for video only (set it to 1 if the video is not needed)
 * 
 *     for(int l = 1; l < frames_per_image; l++)
 *         plot_process += "!cp plot/out" + QString::number(iter) + ".jpeg plot/out" + QString::number(iter + l) + ".jpeg\n";
 * 
 *     fprintf(gp, plot_process.toStdString().c_str());
 *     fflush(gp);
 * 
 *     fprintf(gp, "unset multiplot\n");
 * 
 *     pclose(gp);
 * 
 *     iter += frames_per_image;
 * 
 * }
 * 
@endcode
 <a name="plain-Functionrun_Shoham"></a>
@code
 * template<typename T, typename blas>
 * void step7::ShohamEM<T, blas>::run_Shoham(){
 * 
 *     std::fstream file_like;
 *     file_like.open("output/likelihood.dat", std::ios::out);
 * 
 *     initialize_Shoham();
 * 
 *     dump();
 * 
 *     int iter = 0;
 *     T v_old;
 *     T L_old;
 * 
 *     while(g >= g_min){
 * 
 *         do{
 * 
 *             for(int o = 0; o < g; o++)
 *                 std::cout << pi[o];
 * 
 *             v_old = v;
 *             L_old = L;
 * 
 *             e_step_Shoham();
 * 
 *             m_step_Shoham();
 * 
 *             assemble_P_ij(gamma, v, p, n_data_points, g,
 *                           Delta.array().val(),
 *                           delta_ij.array().val(),
 *                           P_ij.array().val());
 * 
 *             pi_d  = pi;
 * 
 *             Vector c_d;
 * 
 *             c_d = P_ij * pi_d;
 * 
 *             double log_1 = 0.;
 *             for(int j = 0; j< g; j++)
 *                 log_1 += std::log(n_data_points*pi[j]/12);
 * 
 *             double loglike_3 = (N/2)*log_1+(g/2)+(g*(N+1)/2);
 * 
 *             T sum;
 *             sum = sum_array(n_data_points, c_d.array().val(), true);
 * 
 *             L = sum - loglike_3;
 * 
 *             printf("###################loglike = %f ########################\n", L);
 * 
 *         } while(std::fabs(L_old - L) >= 0.1 || std::fabs(v - v_old) >= 0.01);
 * 
 *         dump();
 * 
 *         file_like << L << " "<<L_max <<" "<<v<< " " << std::fabs(L_old - L) << " "<<std::fabs(v - v_old)<<" "<<g<<"\n";
 *         flush(file_like);
 * 
 *         if(L >= L_max + L_max * 0.09 && g > g_min){
 * 
 *             L_max = L;
 * 
 *             pi_opt.resize(pi.size());
 *             std::copy(pi.begin(), pi.end(), pi_opt.begin());
 *             mu_opt = mu;
 *             Sigma_opt.resize(Sigma.size());
 *             std::copy(Sigma.begin(), Sigma.end(), Sigma_opt.begin());
 *             v_opt = v;
 *             u_ij_opt = u_ij;
 *             z_ij_opt = z_ij;
 *             P_ij_opt = P_ij;
 * 
 *             int min = 0;
 * 
 *             for(int i = 1; i < pi.size(); i++)
 *                 if(pi[i] < pi[min])
 *                     min = i;
 * 
 *             pi[min] = 0;
 * 
 *             purge_step_Shoham();
 * 
 *         }
 *         else
 *             break;
 * 
 *         iter++;
 * 
 *     }
 * 
 *     g++;
 *     pi = pi_opt;
 *     mu = mu_opt;
 *     Sigma = Sigma_opt;
 *     v = v_opt;
 *     u_ij = u_ij_opt;
 *     z_ij = z_ij_opt;
 *     P_ij = P_ij_opt;
 * 
 *     dump();
 * 
 *     file_like.close();
 * 
 * }
 * 
@endcode
 <a name="plain-Functioninitialize_Shoham"></a>
@code
 * template<typename T, typename blas>
 * void step7::ShohamEM<T, blas>::initialize_Shoham()
 * {
 *     it = 0;
 *     dealii::IdentityMatrix I(p);
 * 
 *     mu = mean; // cluster centers previous calculated by k-means
 * 
 * 
 *     g_min = 1;
 * 
 *     pi.resize(g);
 *     for(int i = 0; i < g; i++)
 *         pi[i] = 1./g;
 * 
 *     Sigma.resize(g);
 *     for(int i = 0; i < g; i++)
 *         Sigma[i] = I;
 * 
 *     v = 50;
 * 
 *     L_max = -1. * std::numeric_limits<T>::max(); // -infinity
 * 
 *     N = 0.1;//p;
 * 
 *     u_ij.reinit(n_data_points, g);
 *     z_ij.reinit(n_data_points, g);
 *     P_ij.reinit(n_data_points, g);
 * 
 *     delta_ij.reinit(n_data_points, g);
 *     Delta.reinit(g);
 * 
 *     Sigma_inv.resize(g);
 *     for(int i = 0; i < g; i++)
 *         Sigma_inv[i] = I;
 * 
 *     Sigma_dec.resize(g);
 *     for(int i = 0; i < g; i++)
 *         Sigma_dec[i] = I;
 * 
 *     M.reinit(n_data_points, p);
 * 
 *     iter = 0;
 * 
 * 
 * }
 * 
@endcode
 <a name="plain-Functione_step_Shoham"></a>
@code
 * template<typename T, typename blas>
 * void step7::ShohamEM<T, blas>::e_step_Shoham()
 * {
 * 
 *     Matrix tmp(p, n_data_points);
 *     Vector c_d;
 *     std::vector<T> Delta_h(g);
 * 
 *     gamma = tgamma((v+p)/2) / (tgamma(v/2)*(std::pow((M_PI*v),(p/2.))));
 * 
 * 
 *     for(int j=0;j<g;j++)
 *         Delta_h[j] = 1.;
 * 
 *     for(int k = 0; k < g; k++){
 * 
 *         subtract_array_from_matrix(n_data_points, p, k, g,
 *                                    data.array().val(),
 *                                    mu.array().val(),
 *                                    M.array().val());
 * 
 *         tmp = Sigma_inv[k] * tr(M);
 * 
 *         MM_scalar(n_data_points, p,
 *                   M.array().val(),
 *                   tmp.array().val(),
 *                   delta_ij.array().val() + k * n_data_points);
 *     }
 * 
 * 
 *     FullMatrixAccessor delta_ij_h(n_data_points, g, true);
 *     delta_ij_h = delta_ij;
 *     dealii::Vector<T> delta_ij_array(n_data_points * g);
 * 
 * 
 *     for(int i = 0; i <n_data_points; i++){ // transpose
 *         for(int j = 0; j < g; j++){
 *             delta_ij_array(j*n_data_points+i) = delta_ij(i,j);
 *         }
 *     }
 * 
 *     typename dealii::Vector<T>::const_iterator max_delta = std::max_element(delta_ij_array.begin(), delta_ij_array.end());
 *     typename dealii::Vector<T>::const_iterator min_delta = std::min_element(delta_ij_array.begin(), delta_ij_array.end());
 * 
 *     T interval = *max_delta - *min_delta;
 * 
 *     for(int i = 0; i<n_data_points*g; i++){
 *         delta_ij_array(i) -= *min_delta;
 *     }
 *     delta_ij_array /= interval;
 * 
 * 
 *     const unsigned int n_intervals = 50;
 *     dealii::Histogram Histogram;
 *     Histogram.evaluate(delta_ij_array, n_intervals, Histogram.linear);
 * 
 *     std::string delta_gram = "output/delta_gram.dat";
 *     std::ofstream dg(delta_gram.c_str());
 *     Histogram.write_gnuplot(dg);
 * 
 *     Delta = Delta_h;
 * 
 *     printf("gamma = %f\n", gamma);
 * 
 *     assemble_P_ij(gamma, v, p, n_data_points, g,
 *                   Delta.array().val(),
 *                   delta_ij.array().val(),
 *                   P_ij.array().val());
 * 
 *     pi_d  = pi;
 * 
 *     c_d = P_ij * pi_d;
 * 
 *     if(it++==0)
 *         assemble_z_ij(n_data_points, g,
 *                       P_ij.array().val(),
 *                       c_d.array().val(),
 *                       pi_d.array().val(),
 *                       z_ij.array().val());
 * 
 * 
 * 
 * 
 *     assemble_u_ij(v, p, n_data_points, g,
 *                   delta_ij.array().val(),
 *                   u_ij.array().val());
 * 
 * }
 * 
@endcode
 <a name="plain-Functionm_step_Shoham"></a>
@code
 * template<typename T, typename blas>
 * void step7::ShohamEM<T, blas>::m_step_Shoham()
 * {
 * 
 *     double digamma;
 *     Matrix P_ij_tmp = P_ij;
 *     Matrix z_ij_tmp = z_ij;
 *     Matrix q_ij(n_data_points, g);
 *     dealii::IdentityMatrix I(p);
 * 
 *     digamma = boost::math::digamma((p+v)/2.);
 * 
 *     double y = calculate_y(n_data_points,
 *                            g,
 *                            v,
 *                            digamma,
 *                            z_ij.array().val(),
 *                            delta_ij.array().val(),
 *                            u_ij.array().val());
 * 
 *     y /= -n_data_points;
 * 
 *     y = std::max(y, 2.);
 *     printf("y = %f \n", y);
 * 
 *     z_ij_tmp = z_ij;
 *     T max, sum_pi = 2., tmp_sum[g];
 * 
 *     for(int j = 0; j<g; j++){
 *         tmp_sum[j] = 0.;
 *     }
 * 
 *     int rel[g];
 * 
 *     for(int i = 0; i < g; i++) rel[i] = 0;
 * 
 *     for(int i = 0; i < n_data_points; i++){
 * 
 *         int max = 0;
 * 
 *         for(int j = 1; j < g; j++)
 *             if(z_ij(i, j) > z_ij(i, max))
 *                 max = j;
 * 
 *         rel[max]++;
 * 
 *     }
 * 
 *     int s = 0;
 * 
 *     for(int i = 0; i < g; i++)
 *         s+=rel[i];
 * 
 *     for(int i = 0; i < g; i++)
 *         pi[i] = ((double) rel[i]+1) / (double) s;
 * 
 *     for(int i = 0; i < g; i++)
 *         std::cout << pi[i] << ", ";
 * 
 *     purge_step_Shoham();
 * 
 *     std::cout << "\n\n";
 * 
 *     for(int i = 0; i < g; i++)
 *         std::cout << pi[i] << ", ";
 * 
 *     fflush(stdout);
 * 
 *     for(int i = 0; i < g; i++)
 *         if(pi[i] < 0.001) exit(0);
 * 
 *     T zxu[g];
 *     for(int i = 0; i<g; i++){
 *         zxu[i] = 0.;
 *     }
 * 
 *     sum_array2D(n_data_points, g,
 *                 z_ij.array().val(),
 *                 u_ij.array().val(),
 *                 zxu);
 * 
 *     MxM(n_data_points, g,
 *         z_ij.array().val(),
 *         u_ij.array().val(),
 *         q_ij.array().val());
 * 
 *     mu = tr(q_ij) * data;
 * 
 * 
 *     for(int i = 0; i<g; i++){
 *         MatrixSubRow mu_sub(mu,i,0);
 *         mu_sub *= 1./zxu[i];
 *     }
 * 
 * 
 *     M.reinit(n_data_points, p);
 * 
 *     for(int j = 0; j<g; j++){
 *         SMX(n_data_points,
 *             p,
 *             q_ij.array().val(),
 *             zxu[j],
 *             j,
 *             data.array().val(),
 *             mu.array().val(),
 *             M.array().val()
 *             );
 * 
 *         Sigma[j] = tr(M)*M;
 *     }
 * 
 * 
 *     int Matrix_Size = ((p+step1::DEFAULT_TILE_SIZE-1)/step1::DEFAULT_TILE_SIZE)*step1::DEFAULT_TILE_SIZE;
 *     for(int k=0; k < g; k++){
 * 
 * 
 *         Matrix cholMatrix;
 *         dealii::IdentityMatrix eye(Matrix_Size);
 *         cholMatrix = eye;
 * 
 *         SubMatrix Sub_cholMatrix(cholMatrix, 0, p, 0 , p);
 * 
 *         Sub_cholMatrix = Sigma[k];
 * 
 *         step1::Kernels<T> fac_backend;
 * 
 *         fac_backend.cholesky.blackbox(cholMatrix.array().val(), cholMatrix.n_cols(), cholMatrix.leading_dim);
 * 
 *         SubMatrix Sub_Sigma(Sigma_dec[k], 0,0);
 * 
 *         Sub_Sigma = cholMatrix;
 * 
 *         T det = 0.;
 * 
 *         det = Det(p, Matrix_Size, Sub_cholMatrix.array().val());
 *         Delta.add(k, det);
 * 
 *         Matrix Identity = I;
 *         Matrix X(p,p);
 * 
 *         Sigma_inv[k] = Sigma_dec[k];
 *         const SubMatrix l_d(Sigma_dec[k],0,0);
 *         const SubMatrix rhs(Identity, 0,0);
 *         SubMatrix x_d(X,0,0);
 *         SubMatrix inv(Sigma_inv[k],0,0);
 * 
 *         RightMSolveTr<T, blas> L_X(trSub(l_d),x_d);
 *         L_X = rhs;
 * 
 *         RightMSolve<T, blas> L_B(l_d, inv);
 *         L_B = x_d;
 * 
 *     }
 * 
 *     double errfunc = boost::math::erf(0.6594*log(2.1971/(y+log(y)-1)));
 *     v = (2/(y+log(y)-1))+0.0416*(1+errfunc);
 * 
 *     printf("n_dofs_t_distro = %f\n", v);
 * 
 * }
 * 
@endcode
 <a name="plain-Functionpurge_step_Shoham"></a>
@code
 * template<typename T, typename blas>
 * void step7::ShohamEM<T, blas>::purge_step_Shoham()
 * {
 * 
 *     for(int i = 0; i < pi.size(); i++){
 * 
 *         if(pi[i] < 0.0001 && g>g_min){
 * 
 *             pi[i] = pi[g-1];
 *             pi.resize(g-1);
 * 
 *             Sigma[i] = Sigma[g-1];
 *             Sigma.resize(g-1);
 * 
 *             MatrixSubCol subPij_target(P_ij, 0, i);
 *             MatrixSubCol subPij_source(P_ij, 0, g-1);
 *             subPij_target = subPij_source;
 * 
 *             Matrix new_Pij(n_data_points, g-1);
 * 
 *             for(int j = 0; j < g-1; j++){
 * 
 *                 MatrixSubCol to(new_Pij, 0, j);
 *                 MatrixSubCol from(P_ij, 0, j);
 * 
 *                 to = from;
 * 
 *             }
 * 
 *             P_ij = new_Pij;
 * 
 *             MatrixSubCol subDeltaij_target(P_ij, 0, i);
 *             MatrixSubCol subDeltaij_source(P_ij, 0, g-1);
 *             subDeltaij_target = subDeltaij_source;
 * 
 *             Matrix new_Deltaij(n_data_points, g-1);
 * 
 *             for(int j = 0; j < g-1; j++){
 * 
 *                 MatrixSubCol to(new_Deltaij, 0, j);
 *                 MatrixSubCol from(delta_ij, 0, j);
 * 
 *                 to = from;
 * 
 *             }
 * 
 *             delta_ij = new_Deltaij;
 * 
 *             MatrixSubRow subMu_target(mu, i, 0);
 *             MatrixSubRow subMu_source(mu, g-1, 0);
 *             subMu_target = subMu_source;
 * 
 *             Matrix new_mu(g-1, p);
 * 
 *             for(int j = 0; j < g-1; j++){
 * 
 *                 MatrixSubRow to(new_mu, j, 0);
 *                 MatrixSubRow from(mu, j, 0);
 * 
 *                 to = from;
 * 
 *             }
 * 
 *             mu = new_mu;
 * 
 *             MatrixSubCol subUij_target(u_ij, 0, i);
 *             MatrixSubCol subUij_source(u_ij, 0, g-1);
 *             subUij_target = subUij_source;
 * 
 *             Matrix new_Uij(n_data_points, g-1);
 * 
 *             for(int j = 0; j < g-1; j++){
 * 
 *                 MatrixSubCol to(new_Uij, 0, j);
 *                 MatrixSubCol from(u_ij, 0, j);
 * 
 *                 to = from;
 * 
 *             }
 * 
 *             u_ij = new_Uij;
 * 
 *             MatrixSubCol subZij_target(z_ij, 0, i);
 *             MatrixSubCol subZij_source(z_ij, 0, g-1);
 *             subZij_target = subZij_source;
 * 
 *             Matrix new_Zij(n_data_points, g-1);
 * 
 *             for(int j = 0; j < g-1; j++){
 * 
 *                 MatrixSubCol to(new_Zij, 0, j);
 *                 MatrixSubCol from(z_ij, 0, j);
 * 
 *                 to = from;
 * 
 *             }
 * 
 *             z_ij = new_Zij;
 * 
 *             Delta.set(i, Delta(g-1));
 *             Vector newDelta;
 *             newDelta.reinit(g-1);
 * 
 *             for(int j = 0; j < g-1; j++)
 *                 newDelta.set(j, Delta(j));
 * 
 *             Delta = newDelta;
 * 
 *             Sigma_inv[i] = Sigma_inv[g-1];
 *             Sigma_inv.resize(g-1);
 * 
 *             Sigma_dec[i] = Sigma_dec[g-1];
 *             Sigma_dec.resize(g-1);
 * 
 *             g--;
 * 
 *             if(pi[i] == 0) // copied one is zero
 *                 i--;
 * 
 *             purge_step_Shoham();
 * 
 *             break;
 * 
 *         }
 * 
 *     }
 * 
 * }
 * 
@endcode
 <a name="plain-DestructorShohamEM"></a>
@code
 * template<typename T, typename blas>
 * step7::ShohamEM<T, blas>::~ShohamEM()
 * {
 *     BW::Shutdown();
 * }
 * 
@endcode
 <a name="plain-ConstructorKernelTest"></a>
@code
 * template<typename T, typename blas>
 * step7::KernelTest<T, blas>::KernelTest(int n_data_points, int dim, int n_clusters)
 *     :
 *       n_data_points(n_data_points), dim(dim),
 *       n_clusters(n_clusters)
 * {
 *     BW::Init();
 * }
 * 
@endcode
 <a name="plain-Functionrun"></a>
@code
 * template<typename T, typename blas>
 * void step7::KernelTest<T, blas>::run()
 * {
 *     int Matrix_Size = ((dim+step1::DEFAULT_TILE_SIZE-1)/step1::DEFAULT_TILE_SIZE)*step1::DEFAULT_TILE_SIZE;
 * 
 *     FullMatrixAccessor Test_Matrix(n_data_points, dim, true);
 *     FullMatrixAccessor Test_Matrix2(n_clusters, dim, true);
 *     FullMatrixAccessor Test_Matrix3(n_data_points, dim, true);
 *     FullMatrixAccessor Test_Matrix4(dim, n_data_points, true);
 * 
 *     FullMatrixAccessor Test_Matrix_ij(n_data_points, n_clusters, true);
 *     FullMatrixAccessor Test_Matrix_ij2(n_data_points, n_clusters, true);
 *     FullMatrixAccessor Test_Matrix_ij3(n_data_points, n_clusters, true);
 * 
 *     std::vector<T> Delta_h(n_clusters);
 *     std::vector<T> c_h(n_data_points);
 *     std::vector<T> pi_h(n_clusters);
 *     Vector Delta, c, pi;
 * 
 *     T zxu[n_clusters];
 *     for(int i = 0; i<n_clusters; i++){
 *         zxu[i] = 0.;
 *     }
 * 
 *     for(int i = 0; i< n_data_points; i++){
 *         for(int j=0; j<dim; j++){
 *             Test_Matrix(i,j) = (i+1)/10.;
 *             Test_Matrix3(i,j) = 0.;
 *             Delta_h[j] = 1.;
 *             c_h[i] = 1.;
 *             pi_h[j] = 1/n_clusters;
 *         }
 *     }
 * 
 *     for(int i = 0; i< n_clusters; i++){
 *         for(int j=0; j<dim; j++){
 *             Test_Matrix2(i,j) = (i+1)/10.;
 *         }
 *     }
 * 
 *     for(int i = 0; i< dim; i++){
 *         for(int j=0; j<n_data_points; j++){
 *             Test_Matrix4(i,j) = (i+1)/10.;
 *         }
 *     }
 * 
 *     for(int i = 0; i< n_data_points; i++){
 *         for(int j=0; j<n_clusters; j++){
 *             Test_Matrix_ij(i,j) = (i+1)/10.;
 *             Test_Matrix_ij2(i,j) = (i+1)/10.;
 *             Test_Matrix_ij3(i,j) = (i+1)/10.;
 *         }
 *     }
 * 
 *     Matrix Test, Test2, Test3, Test4, Test_ij, Test_ij2, Test_ij3;
 *     Test = Test_Matrix;
 *     Test2 = Test_Matrix2;
 *     Test3 = Test_Matrix3;
 *     Test4 = Test_Matrix4;
 *     Test_ij = Test_Matrix_ij;
 *     Test_ij2 = Test_Matrix_ij2;
 *     Test_ij3 = Test_Matrix_ij3;
 *     Delta = Delta_h;
 *     c = c_h;
 *     pi = pi_h;
 * 
 *     subtract_array_from_matrix(n_data_points, dim, 1, n_clusters,
 *                                Test.array().val(),
 *                                Test2.array().val(),
 *                                Test3.array().val());
 * 
 *     printf("dim = %d\n", dim);
 * 
 * 
 *     MM_scalar(n_data_points, dim,
 *               Test3.array().val(),
 *               Test4.array().val(),
 *               Test.array().val());
 * 
 *     Det(dim, Matrix_Size, Test3.array().val());
 * 
 *     assemble_P_ij(0.5, 50, dim, n_data_points, n_clusters,
 *                   Delta.array().val(),
 *                   Test_ij.array().val(),
 *                   Test_ij2.array().val());
 * 
 *     assemble_u_ij(50, dim, n_data_points, n_clusters,
 *                   Test_ij.array().val(),
 *                   Test_ij2.array().val());
 * 
 *     assemble_z_ij(n_data_points, n_clusters,
 *                   Test_ij.array().val(),
 *                   c.array().val(),
 *                   pi.array().val(),
 *                   Test_ij2.array().val());
 * 
 *     calculate_y(n_data_points, n_clusters, 0.5, 13,
 *                 Test_ij.array().val(),
 *                 Test_ij2.array().val(),
 *                 Test_ij3.array().val());
 * 
 *     sum_array(n_data_points, c.array().val(), false);
 * 
 *     sum_array2D(n_data_points, n_clusters,
 *                 Test_ij.array().val(),
 *                 Test_ij2.array().val(),
 *                 zxu);
 * 
 *     MxM(n_data_points, n_clusters,
 *         Test_ij.array().val(),
 *         Test_ij2.array().val(),
 *         Test_ij3.array().val());
 * 
 *     SMX(n_data_points, dim,
 *         Test_ij.array().val(),
 *         zxu[2],
 *         2,
 *         Test.array().val(),
 *         Test2.array().val(),
 *         Test3.array().val());
 * }
 * 
@endcode
 <a name="plain-DestructorKernelTest"></a>
@code
 * template<typename T, typename blas>
 * step7::KernelTest<T, blas>::~KernelTest()
 * {
 *     BW::Shutdown();
 * }
 * #endif // CUDA_DRIVER_STEP_7_HH
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #include <iostream>
 * #include <vector>
 * #include <fstream>
 * 
 * #include <QString>
 * 
 * #include <cuda_driver_step-7.h>
 * 
 * #include <deal.II/lac/matrix_out.h>
 * 
 * #include <cuda_driver_step-7.hh>
 * 
 * #include <waveforms_generator.h>
 * 
 * #include <step-4.hh>
 * 
 * #include <cuda_driver_step-5.h>
 * #include <cuda_driver_step-5.hh>
 * 
 * namespace step7 {
 * 
@endcode
 <a name="plain-ClassSimParams"></a>
@code
 * template<typename T>
 * struct SimParams : public step4::QRTestUIParams
 *         #ifdef USE_KMEANS
 *         , step6::SimParam
 *         #endif
 * 
 * {
 * 
 *     typedef step4::QRTestUIParams Base1;
 * 
 * #ifdef USE_KMEANS
 *     typedef step6::SimParam Base2;
 * #endif
 * 
 *     SimParams() : Base1()
 *   #ifdef USE_KMEANS
 *       , Base2()
 *   #endif
 *     {}
 * 
 *     static void declare(dealii::ParameterHandler & prm);
 * 
 *     void get(dealii::ParameterHandler & prm);
 * 
 *     dealii::FullMatrix<T> waveforms;
 *     dealii::FullMatrix<T> check_waves;
 * 
 *     int n_forms, n_rows, n_cols;
 *     const static int forms = 3;
 *     T noise;
 *     std::vector<double> waves_parameters;
 *     std::vector<int> n_waves;
 *     QString number_of_waves;
 *     std::string output_folder;
 * 
 *     unsigned int n_components, max_iter;
 * 
 *     double ev_tol;
 * 
 * public:
 * 
 *     void generate_waveforms();
 * 
 * };
 * 
 * 
@endcode
 <a name="plain-ClassPCA"></a>
@code
 * template <typename T, typename blas> class PCA {
 * 
 * public:
 * 
 *     typedef typename blas_pp<T, blas>::FullMatrixAccessor FullMatrixAccessor;
 *     typedef typename blas_pp<T, blas>::Matrix       Matrix;
 *     typedef          transpose<Matrix>              tr;
 * 
 *     typedef typename blas_pp<T, blas>::SubMatrix    SubMatrix;
 *     typedef typename blas_pp<T, blas>::MatrixSubCol MatrixSubCol;
 *     typedef typename blas_pp<T, blas>::Vector       Vector;
 * 
 * 
 *     PCA(dealii::ParameterHandler & prm);
 *     void run();
 * 
 * 
 *     Matrix Q_x_scores;
 *     Matrix d_scores;
 *     Matrix d_loadings;
 * 
 * private:
 *     void factorize(dealii::FullMatrix<T> &A);
 * 
 *     void check_results(const dealii::FullMatrix<T> & A,
 *                        const step4::CudaQRDecomposition<T, blas>& QRf,
 *                        T elapsed_time);
 * 
 *     void save_results();
 * 
 *     dealii::FullMatrix<T> Q, B, P_t, H;
 * 
 *     SimParams<T> params;
 * 
 *     dealii::ConvergenceTable results_table;
 * 
 *     unsigned int n_successful_measurements;
 * };
 * 
@endcode
 <a name="plain-ClassMyFancySimulation"></a>
@code
 * template<typename T>
 * class MyFancySimulation {
 * 
 * public:
 * 
 *     typedef typename blas_pp<T, blas>::FullMatrixAccessor FullMatrixAccessor;
 *     typedef typename blas_pp<T, blas>::Matrix       Matrix;
 *     typedef          transpose<Matrix>              tr;
 * 
 *     typedef typename blas_pp<T, blas>::SubMatrix    SubMatrix;
 *     typedef typename blas_pp<T, blas>::MatrixSubCol MatrixSubCol;
 *     typedef typename blas_pp<T, blas>::Vector       Vector;
 * 
 *     MyFancySimulation(std::string prm_filename);
 * 
 *     void run();
 * 
 *     void run_pca(dealii::ParameterHandler &prm_handler);
 * 
 *     void run_cluster_analysis();
 * 
 * private:
 * 
 *     SimParams<T> params;
 * 
 *     FullMatrixAccessor data;
 *     FullMatrixAccessor loadings;
 *     FullMatrixAccessor  scores;
 * };
 * 
 * 
@endcode
 <a name="plain-ConstructorPCA"></a>
@code
 * template <typename T, typename blas>
 * step7::PCA<T, blas>::PCA(dealii::ParameterHandler &prm)
 *     :
 *       n_successful_measurements(0)
 * {
 *     params.get(prm);
 * }
 * 
@endcode
 <a name="plain-Functionrun"></a>
@code
 * template <typename T, typename blas>
 * void step7::PCA<T, blas>::run() {
 * 
 *     params.generate_waveforms();
 * 
 *     std::string wforms = "output/wforms.txt";
 *     std::ofstream wf(wforms.c_str());
 *     params.waveforms.print(wf, 12 / *width; avoids that two numbers appear as one in the output* /);
 * 
 *     factorize(params.waveforms);
 * 
 * }
 * 
@endcode
 <a name="plain-Functionfactorize"></a>
@code
 * template <typename T, typename blas>
 * void
 * step7::PCA<T, blas>::factorize(dealii::FullMatrix<T> &A)
 * {
 * 
 *     std::cout  << std::endl << " ---------- FACTORIZE ----------" << std::endl;
 * 
 *     step4::CudaQRDecomposition<T, blas> QR;
 * 
 *     QR.householder(A);
 * 
 *     std::cout  << std::endl << " ---------- DONE ----------" << std::endl;
 * 
 *     FullMatrixAccessor R_tmp;
 *     R_tmp = QR.R();
 * 
 * 
 *     dealii::FullMatrix<T> R_tmp_h(R_tmp.n_rows(), R_tmp.n_cols());
 * 
 *     std::ofstream R_out("output/R_mea.dat");
 * 
 *     std::cout<< "R dimensions = " << R_tmp.n_rows() << " X " << R_tmp.n_cols() << std::endl;
 *     for (unsigned int r = 0; r < R_tmp_h.n_rows(); ++r)
 *     {
 *         for(unsigned int c = 0; c < R_tmp_h.n_cols(); ++c)
 *         {
 *             R_tmp_h(r, c) = R_tmp(r, c);
 *             R_out << R_tmp_h(r, c) << "\t";
 *         }
 *         R_out  << std::endl;
 *     }
 *     std::cout  << std::endl;
 * 
 *     std::cout  << std::endl << " ---------- NIPALS ----------" << std::endl;
 * 
 *     step5::CudaNIPALS<T, blas> nipals(params.n_components,
 *                                       -1/ *params.max_iter* /,
 *                                       params.ev_tol);
 * 
 *     nipals.get_pca(R_tmp_h);
 * 
 *     std::cout  << std::endl << " ---------- DONE ----------" << std::endl;
 * 
 *     nipals.lambda.print(std::cout);
 * 
 * 
 *     Q_x_scores.reinit(params.n_rows,params.n_components);
 * 
 *     SubMatrix sub_Q(const_cast<Matrix &>(QR.Q()),0,params.n_rows,0,params.n_cols);
 * 
 *     SubMatrix scores(const_cast<Matrix &>(nipals.scores()),0,0);
 * 
 *     SubMatrix Q_x_s(const_cast<Matrix &>(Q_x_scores),0,0);
 * 
 * 
 *     d_scores.reinit(nipals.scores().n_rows(), nipals.scores().n_cols());
 *     d_scores = nipals.scores();
 * 
 *     d_loadings.reinit(nipals.loads().n_rows(),nipals.loads().n_cols());
 *     d_loadings = nipals.loads();
 * 
 *     Matrix Q_tmp, S_tmp, P_tmp;
 *     Q_tmp = QR.Q();
 *     S_tmp = nipals.scores();
 *     P_tmp = Q_tmp * S_tmp;
 * 
 *     Q_x_scores = P_tmp;
 * 
 *     std::ofstream PCA_out("output/pca_out.txt");
 *     for(unsigned int r = 0; r < P_tmp.n_rows(); r++){
 *         for(unsigned int c = 0; c < P_tmp.n_cols(); c++)
 *             PCA_out<<P_tmp(r,c)<<"\t";
 *         PCA_out<<"\n";
 *     }
 * 
 * }
 * 
@endcode
 <a name="plain-Functiondeclare"></a>
@code
 * template<typename T>
 * void step7::SimParams<T>::declare(dealii::ParameterHandler &prm)
 * {
 *     Base1::declare(prm);
 * #ifdef USE_KMEANS
 *     Base2::declare(prm);
 * #endif
 * 
 *     prm.enter_subsection("Waveforms parameters");
 * 
 *     prm.declare_entry("Number of waves", "300,300,300",
 *                       dealii::Patterns::Anything(),
 *                       "Number of Waves");
 * 
 *     prm.declare_entry("Noise intensity", "0.2",
 *                       dealii::Patterns::Double(0.,0.9),
 *                       "Noise intensity");
 * 
 *     prm.declare_entry("Number of columns", "128",
 *                       dealii::Patterns::Integer(),
 *                       "Number of columns");
 * 
 *     prm.declare_entry("Center of Wave 1", "0.25",
 *                       dealii::Patterns::Double(0.,0.9),
 *                       "Center of Wave 1");
 * 
 *     prm.declare_entry("Sigma of Wave 1", "0.09",
 *                       dealii::Patterns::Double(0.01,0.09),
 *                       "Sigma of Wave 1 (standard deviation of gaussian)");
 * 
 *     prm.declare_entry("Center of Wave 2", "0.85",
 *                       dealii::Patterns::Double(0., 1.),
 *                       "Center of Wave 2");
 * 
 *     prm.declare_entry("Center 1 of Wave 3", "0.6",
 *                       dealii::Patterns::Double(0.,0.9),
 *                       "Center 1 of Wave 3");
 * 
 *     prm.declare_entry("Center 2 of Wave 3", "0.65",
 *                       dealii::Patterns::Double(0.,0.9),
 *                       "Center 2 of Wave 3");
 * 
 *     prm.declare_entry("Sigma 1 of Wave 3", "0.05",
 *                       dealii::Patterns::Double(0.01,0.09),
 *                       "Sigma 1 of Wave 3 (Lorenz)");
 * 
 *     prm.declare_entry("Sigma 2 of Wave 3", "0.06",
 *                       dealii::Patterns::Double(0.01,0.09),
 *                       "Sigma 2 of Wave 3");
 * 
 *     prm.declare_entry ("Output folder",
 *                        "output",
 *                        dealii::Patterns::Anything(),
 *                        "output folder location");
 * 
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("NIPALS parameters");
 * 
 *     prm.declare_entry("Number of components", "2",
 *                       dealii::Patterns::Integer(),
 *                       "Number of components");
 * 
 *     prm.declare_entry("Maximal number of iterations", "1000",
 *                       dealii::Patterns::Integer(),
 *                       "Maximal number of iterations. "
 *                       "If -1 is given termination of NIPALS iterations "
 *                       "solely depends on the given tolerance");
 * 
 *     prm.declare_entry("Tolerance for eigenvalues",
 *                       "1e-5",
 *                       dealii::Patterns::Double(),
 *                       "Once the relative change of an eigenvalue in the inner "
 *                       "power iteration is les than the given value iteration "
 *                       "stops and the next principal component is computed.");
 * 
 *     prm.leave_subsection();
 * 
 * }
 * 
@endcode
 <a name="plain-Functionget"></a>
@code
 * template<typename T>
 * void step7::SimParams<T>::get(dealii::ParameterHandler &prm)
 * {
 * 
 *     this->Base1::get(prm);
 * 
 * #ifdef USE_KMEANS
 *     this->Base2::get(prm);
 * #endif
 * 
 *     prm.enter_subsection("Waveforms parameters");
 * 
 *     number_of_waves = QString(prm.get("Number of waves").c_str());
 * 
 *     noise = prm.get_double("Noise intensity");
 * 
 *     n_cols = prm.get_integer("Number of columns");
 * 
 *     waves_parameters.push_back(prm.get_double("Center of Wave 1"));
 * 
 *     waves_parameters.push_back(prm.get_double("Sigma of Wave 1"));
 * 
 *     waves_parameters.push_back(prm.get_double("Center of Wave 2"));
 * 
 *     waves_parameters.push_back(prm.get_double("Center 1 of Wave 3"));
 * 
 *     waves_parameters.push_back(prm.get_double("Center 2 of Wave 3"));
 * 
 *     waves_parameters.push_back(prm.get_double("Sigma 1 of Wave 3"));
 * 
 *     waves_parameters.push_back(prm.get_double("Sigma 2 of Wave 3"));
 * 
 *     output_folder = prm.get ("Output folder");
 * 
 *     QDir dir("plot");
 * 
 *     if(!dir.exists()){
 *         bool mkdir = dir.mkpath("./");
 *         AssertThrow(mkdir, dealii::ExcMessage("creating a new Folder failed"));
 *     }
 * 
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("NIPALS parameters");
 * 
 *     this->n_components = prm.get_integer("Number of components");
 * 
 *     this->max_iter = prm.get_integer("Maximal number of iterations");
 * 
 *     this->ev_tol   = prm.get_double("Tolerance for eigenvalues");
 * 
 *     prm.leave_subsection();
 * 
 * }
 * 
@endcode
 <a name="plain-Functiongenerate_waveforms"></a>
@code
 * template<typename T>
 * void step7::SimParams<T>::generate_waveforms()
 * {
 *     QStringList number_list = number_of_waves.split(",", QString::SkipEmptyParts);
 * 
 *     if(number_list.size()<forms)
 *     {
 *         std::cout<<"more waves numbers needed!"<<std::endl;
 *         exit(1);
 *     }
 * 
 *     QStringList::iterator e = number_list.begin();
 *     QStringList::iterator end = number_list.end();
 * 
 *     int i=0, sum = 0;
 *     n_waves.resize(forms);
 *     n_forms = 0;
 * 
 *     for(; e!= end; ++e, ++i)
 *     {
 *         if(i==forms) break;
 *         n_waves[i] = e->toInt();
 *     }
 * 
 *     for(int i = 0; i<forms;i++){
 *         sum += n_waves[i];
 *         if(n_waves[i]!=0) n_forms++;
 *     }
 * 
 *     n_rows = sum;
 * 
 *     waveforms.reinit(n_rows,n_cols);
 *     check_waves.reinit(n_forms,n_cols);
 *     WaveForms::generate_waveforms(n_rows, n_cols, waveforms, n_waves, noise, waves_parameters);
 *     WaveForms::normalize_waveforms(n_rows, n_cols, waveforms);
 *     WaveForms::original_waveforms(n_cols,n_forms,check_waves, waves_parameters, n_waves);
 * 
 *     dealii::MatrixOut matrix_out;
 *     std::string wave_forms = output_folder+"/waveforms.txt";
 *     std::ofstream waf(wave_forms.c_str());
 *     matrix_out.build_patches(waveforms, "waveforms");
 *     matrix_out.write_gnuplot(waf);
 * 
 *     dealii::MatrixOut original_out;
 *     std::string original_forms = output_folder+"/originalforms.txt";
 *     std::ofstream oaf(original_forms.c_str());
 *     original_out.build_patches(check_waves, "originalforms");
 * 
 * }
 * 
 * }
 * 
@endcode
 <a name="plain-ConstructorMyFancySimulation"></a>
@code
 * template<typename T>
 * step7::MyFancySimulation<T>::MyFancySimulation(std::string prm_filename)
 * {
 *     dealii::ParameterHandler param_handler;
 *     SimParams<T>::declare(param_handler);
 * 
 *     param_handler.read_input(prm_filename);
 *     params.get(param_handler);
 * }
 * 
 * 
@endcode
 <a name="plain-Funktionrun"></a>
@code
 * template<typename T>
 * void step7::MyFancySimulation<T>::run()
 * {
 * 
 * 
 * }
 * 
 * template<typename T>
 * void step7::MyFancySimulation<T>::run_pca(dealii::ParameterHandler &prm_handler)
 * {
 * 
 *     std::cout << "Householder-QR mit cublas :\n"
 *               << "------------------------------------------------" << std::endl;
 * 
 *     PCA<T, blas> cublas_pca(prm_handler);
 * 
 *     cublas_pca.run();
 *     data = cublas_pca.Q_x_scores;
 *     loadings = cublas_pca.d_loadings;
 *     scores = cublas_pca.d_scores;
 *     std::cout << "\nHouseholder-QR mit cublas DONE \n"
 *               << "------------------------------------------------\n" << std::endl;
 * 
 * }
 * 
@endcode
 <a name="plain-Funktionrun_cluster_analysis"></a>
@code
 * template<typename T>
 * void step7::MyFancySimulation<T>::run_cluster_analysis()
 * {
 * 
 *     std::cout  << std::endl << " ---------- START SIMULATION ----------" << std::endl;
 * 
 *     std::cout  << std::endl << " ---------- GENERATE WAVEFORMS ----------" << std::endl;
 * 
 *     params.generate_waveforms();
 * 
 *     std::cout  << std::endl << " ---------- DONE ----------" << std::endl;
 * 
 *     std::cout  << std::endl;
 * 
 *     FullMatrixAccessor X(params.waveforms,true);
 * 
 *     FullMatrixAccessor Y(params.check_waves, true);
 * 
 *     dealii::FullMatrix<T> data_matrix2(params.n_rows, params.n_cols);
 *     dealii::FullMatrix<T> data_matrix2t(params.n_cols, params.n_rows);
 *     X.push_to(data_matrix2);
 * 
 *     for(int i = 0; i < params.n_rows; i++){
 *         for(int j = 0; j < params.n_cols; j++){
 *             data_matrix2t(j, i) = data_matrix2(i, j);
 *         }
 *     }
 * 
 *     std::string sdata2 = "output/data22.txt";
 *     std::ofstream d2(sdata2.c_str());
 *     data_matrix2.print(d2, 12 / *width; avoids that two numbers appear as one in the output* /);
 *     std::string sdata2t = "output/gnuplot_inputwaves.txt";
 *     std::ofstream d2t(sdata2t.c_str());
 *     data_matrix2t.print(d2t, 12 / *width; avoids that two numbers appear as one in the output* /);
 * 
 * 
 *     std::cout  << std::endl << " ---------- DONE ----------" << std::endl;
 * 
 * #ifdef USE_KMEANS
 * 
 *     std::cout  << std::endl << " ---------- KMEANS ----------" << std::endl;
 * 
 *     std::cout << "n_clusters= " << params.n_clusters <<std::endl;
 * 
 *     step6::Kmeans<double, cublas> cluster_search(data,params.mmeans,params.max_iter);
 *     cluster_search.initial_computation(params.method);
 *     cluster_search.iterate(params.method,params.max_iter);
 * 
 *     std::string smmeans = "output/mmeans.txt";
 *     std::ofstream j(smmeans.c_str());
 *     params.mmeans.print(j, 12 / *width; avoids that two numbers appear as one in the output* /);
 * 
 *     std::cout  << std::endl << " ---------- DONE ----------" << std::endl;
 * 
 * #endif
 * 
 *     dealii::ParameterHandler prm_debug_handler;
 * 
 *     FullMatrixAccessor means(params.mmeans, true);
 * 
 *     std::cout  << std::endl << " ---------- SHOHAM ----------" << std::endl;
 * 
 *     step7::ShohamEM<T, cublas> Test(data, means); // data, means
 * 
 *     Test.run_Shoham();
 * 
 *     std::cout  << std::endl << " ---------- DONE WITH SHOHAM----------" << std::endl;
 * 
 * #ifdef USE_KMEANS
 *     std::string smeans = "output/means.txt";
 *     std::ofstream m(smeans.c_str());
 *     params.mmeans.print(m, 12 / *width; avoids that two numbers appear as one in the output* /);
 * #endif
 * 
 *     std::string oforms = "output/oforms.txt";
 *     std::ofstream of(oforms.c_str());
 *     params.check_waves.print(of, 12 / *width; avoids that two numbers appear as one in the output* /);
 * 
 * #ifdef USE_KMEANS
 *     dealii::MatrixOut means_out;
 *     std::string means_forms = params.output_folder+"/means-waveforms.txt";
 *     std::ofstream mf(means_forms.c_str());
 *     means_out.build_patches(params.mmeans, "waveforms");
 *     means_out.write_gnuplot(mf);
 * #endif
 * 
 * }
 * 
@endcode
 <a name="plain-Funktionrun_simulation"></a>
@code
 * void run_simulation(int / *argc* /, char *argv[])
 * {
 *     dealii::ParameterHandler prm_handler;
 * 
 * #ifdef USE_SINGLE_PRECISION
 *     typedef float NumberType;
 * #else
 *     typedef double NumberType;
 * #endif
 * 
 *     typedef step7::SimParams<NumberType> SimParams;
 * 
 *     SimParams::declare(prm_handler);
 * 
 *     std::string prm_filename(argv[0]);
 *     prm_filename += ".prm";
 *     prm_handler.read_input (prm_filename);
 * 
 *     prm_filename += ".log";
 *     std::ofstream log_out_text(prm_filename.c_str());
 *     prm_handler.print_parameters (log_out_text,
 *                                   dealii::ParameterHandler::Text);
 * 
 *     SimParams params;
 * 
 *     params.get(prm_handler);
 *     int DevNo = 0;
 * 
 *     cudaSetDevice(DevNo);           // select CUDA device
 * 
 *     global_data.current_device_id = DevNo;
 * 
 *     global_data.cublanc_gpu_info(); // obtain informations about used GPU
 * 
 * #ifdef USE_SINGLE_PRECISION
 *     std::cout << "Simulation with float:\n"
 *               << "------------------------------------------------" << std::endl;
 * 
 *     step7::MyFancySimulation<float> simu_float(prm_filename);
 * 
 *     simu_float.run_pca(prm_handler);
 * 
 *     simu_float.run_cluster_analysis();
 * 
 *     std::cout << "\nSimulation with float DONE \n"
 *               << "------------------------------------------------\n" << std::endl;
 * #else
 *     std::cout << "Simulation with double:\n"
 *               << "------------------------------------------------" << std::endl;
 * 
 *     fflush(stdout);
 * 
 *     step7::MyFancySimulation<double> simu_double(prm_filename);
 * 
 *     simu_double.run_pca(prm_handler);
 * 
 *     simu_double.run_cluster_analysis();
 * 
 *     std::cout << "\nSimulation with double DONE \n"
 *               << "------------------------------------------------\n" << std::endl;
 * #endif
 * 
 * }
 * 
 * 
@endcode
 <a name="plain-Funktionrun_kerneltest"></a>
@code
 * void run_kerneltest(int / *argc* /, char *argv[])
 * {
 * 
 *     int DevNo = 0;
 * 
 *     cudaSetDevice(DevNo);
 * #ifdef USE_SINGLE_PRECISION
 *     typedef float NumberType;
 * #else
 *     typedef double NumberType;
 * #endif
 * 
 *     static const int n_data_points = 100000;
 *     static const int n_sampling_points = 64;
 *     static const int n_clusters = 10;
 * 
 *     step7::KernelTest<NumberType, cublas> Test(n_data_points,
 *                                                n_sampling_points,
 *                                                n_clusters);
 * 
 *     Test.run();
 * 
 * }
 * 
 * namespace step7 {
 * 
 * }// namespace step7 END
 * 
@endcode
 <a name="plain-Funktionmain"></a>
@code
 * int main(int argc, char *argv[])
 * {
 *     cudaGetDeviceCount(&global_data.n_CUDA_devices);
 *     std::cout
 *             << "N available CUDA devices: "
 *             << global_data.n_CUDA_devices << std::endl;
 * 
 *     run_simulation(argc, argv);
 * 
 * }
 * 
 * 
 * 
 @endcode
 */
