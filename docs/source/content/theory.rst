A coupled multi-domain method for 1D hemodynamics
#################################################

Flow and pressure waves, originating due to the contraction of the heart, propagate along the deformable vessels and reflect due to tapering, branching, and other discontinuities. The size and complexity of the cardiovascular system necessitate a multidomain approach, with “upstream” regions of interest (large arteries) coupled to reduced-order models of “downstream” vessels. Previous efforts to couple upstream and downstream domains have included specifying resistance and impedance outflow boundary conditions for the nonlinear one-dimensional wave propagation equations and iterative coupling between three-dimensional and one-dimensional numerical methods. We have developed a new approach to solve the one-dimensional nonlinear equations of blood flow in elastic vessels utilizing a space-time finite element method with GLSstabilization for the upstream domain, and a boundary term to couple to the downstream domain. 

The outflow boundary conditions are derived following an approach analogous to the Dirichlet-to-Neumann (DtN) method. 

In the downstream domain, we solve simplified 0D/1D equations to derive relationships between pressure and flow accommodating periodic and transient phenomena with a consistent formulation for different boundary condition types. In this chapter, we also present a new boundary condition that accommodates transient phenomena based on a Green’s function solution of the linear, damped wave equation in the downstream domain. 

The mathematical formulation, the numerical derivation and results are presented in the next sections.

We present here the different steps that are required to develop the coupled multidomain method from the strong form in the original domain to the variational form in the computational domain that includes the information of the analytical domain.

Strong form
-----------

The one-dimensional equations for the flow of a Newtonian, incompressible fluid in a deforming, elastic domain consist of the continuity equation, a single axial momentum balance equation, a constitutive equation, and suitable initial and boundary conditions. The governing equations are derived in a general form by Hughes [2]_ and Hughes and Lubliner [1]_. The partial differential equations for mass and momentum balance are given by (z is the axial coordinate):

.. math:: 

   \frac{\partial S}{\partial t} + \frac{\partial Q}{\partial z} = -\psi
          
   \frac{\partial Q}{\partial t} + \frac{\partial}{\partial z}\left[(1 + \delta)\,\frac{Q^2}{S}\right] + \frac{S}{\rho}\,\frac{\partial p}{\partial z} = S\,f + N\,\frac{Q}{S} + \nu\frac{\partial^2 Q}{\partial z^2}

The primary variables are the cross-sectional area :math:`S`, the pressure :math:`p`, and the volumetric flow rate :math:`Q`. The density of the fluid is given by ρ (assumed constant), the external force by :math:`f`, the kinematic viscosity by :math:`\nu` (assumed constant) and :math:`\psi` is an outflow function (taken to be zero for impermeable vessels). 
The variables :math:`\delta` and :math:`N` are defined by the choice of a profile function for the velocity over the cross-section. Here we choose a time-varying, parabolic flow profile, thus [1]_:

.. math::
   :label: eq3

   \delta = \frac{1}{3},\quad N = -8\,\pi\,\nu

The governing equations are of mixed parabolic-hyperbolic type, but have mainly a hyperbolic nature since the diffusive term is small. We thus impose one boundary condition at each inlet/outlet by specifying values of the primary variables or a relationship between them. 

The flow rate is typically specified at the inlet(s) (:math:`\Gamma_{in}`), but the inlet boundaries can accommodate the same types of boundary conditions as will be subsequently discussed for the outlets:

.. math::

   Q(z,t) = Q_{in}(t),\quad z\in\Gamma_{in}

The initial conditions for this problem are given by (where :math:`S^0(z)`, :math:`Q^0(z)` and :math:`p^0(z)` are prescribed functions):

.. math::

   S(z,0) = S^0(z),Q(z,0) = Q^0(z)\,\text{and}\, p(z,0) = p^0(z)

In order to complete the above system, we need to introduce a constitutive relationship. An elastic model is assumed, which relates the pressure to the cross-sectional area as follows:

.. math::

   p(z,t) = \tilde{p}[S(z,t),z,t]

and its inverse function

.. math::
   :label: eq7

   S(z,t) = \tilde{S}[p(z,t),z,t]

The particular constitutive relationship that we have used is that proposed by Olufsen [3]_:

.. math::
   :label: eq8

   \tilde{p}(S,z) = p^0(z) + \frac{4}{3}\,\frac{E\,h}{r^0(z)}\,\left(1 - \sqrt{\frac{S^0(z)}{S(z,t)}}\right)

here the Young’s modulus :math:`E` and the wall thickness :math:`h` relate to the radius :math:`r^0 = \sqrt{S^0(z)/\pi}`:

.. math::
   
   \frac{E\,h}{r^0(z)} = k_1\,\exp{k_2\,r^0(z)} + k_3

In this relationship, :math:`k_1`, :math:`k_2`, and :math:`k_3` are empirically-derived constants with values in CGS units of 2x107 g⋅ s-2⋅ cm-1, -22.53 cm-1, and 8.65x105 g⋅ s-2⋅ cm-1, respectively. 
Here we use a constant initial pressure :math:`p^0(z) = p^0`. By noting that the pressure gradient can be expanded as

.. math::

  \frac{\partial p}{\partial z} = \frac{\partial\tilde{p}}{\partial S}\,\frac{\partial S}{\partial z} + \frac{\partial\tilde{p}}{\partial z}

we can rewrite the system of partial differential equations in the following quasi-linear conservative form:

.. math::
   
   \frac{\partial\mathbf{U}}{\partial t} + \frac{\partial\mathbf{F}}{\partial z} - \mathbf{K}\,\frac{\partial^2\mathbf{U}}{\partial z^2} = G,\,\text{or}\quad\frac{\partial\mathbf{U}}{\partial t} + \frac{\partial\mathbf{F}}{\partial z} - \mathbf{K}\,\frac{\partial^2\mathbf{U}}{\partial z^2} = \mathbf{C}_F\,\mathbf{U}

where

.. math::
   :label: eq13

   \mathbf{U} = 
   \begin{bmatrix}
   U_1\\
   U_2
   \end{bmatrix}
    = 
   \begin{bmatrix}
   S\\
   Q
   \end{bmatrix}

   \mathbf{F} = 
   \begin{bmatrix}
   U_2\\
   (1 + \delta)\,\frac{U_2^2}{U_1} + \frac{1}{\rho}\,\int_{p^0}^{p(z,t)}\tilde{S}(p,z,t)\,dp
   \end{bmatrix},
   \quad
   \mathbf{K} = 
   \begin{bmatrix}
   0 & 0\\
   0 & \nu\\
   \end{bmatrix}.

   \mathbf{G} = 
   \begin{bmatrix}
   -\psi\\
   U_1\,f + N\,\frac{U_2}{U_1} + \int_{p^0}^{p}\frac{1}{\rho}\,\frac{\partial\tilde{S}(p,z,t)}{\partial z}\,dp
   \end{bmatrix},\quad
   \mathbf{C}_F = 
   \begin{bmatrix}
   -\frac{\psi}{U_1} & 0\\
   f + \frac{1}{U_1}\,\int_{p^0}^{p}\frac{1}{\rho}\,\frac{\partial\tilde{S}(p,z,t)}{\partial z}\,dp & \frac{N}{U_1}
   \end{bmatrix}.
   

The motivation to work with the conservative form rather than the advective form as in previous work [4]_, is to be able to integrate by part the convective term and obtain a flux (a boundary integral) through which the multidomain coupling can be performed.

Note that in the advective form, the only term that can easily be integrated by parts is the longitudinal viscous term, which is very small and often neglected in one-dimensional theory. Thus, the main difference between the two forms is the treatment of the boundary conditions. 

In the present conservative formulation, boundary conditions are prescribed in a natural way. In contrast, in the advective form, boundary conditions are enforced in an essential way: the equation for the corresponding dof is replaced by an equation representing the boundary condition.

Weak form
---------

The weak formulation of the initial boundary value problem is given as follows with :math:`\Omega = [0, L]` : find :math:`\mathbf{U}` in :math:`\mathcal{V} = \left\{\mathbf{U}:\Omega\times (0,T)\rightarrow\mathbb{R}^2\,|\,\mathbf{U}(z,t)\in H_0^1\right\}` such that :math:`\forall\,\mathbf{W} = \left[W_1\,W_2\right]^T\in\mathcal{V}`,

.. math::

   \begin{eqnarray}
   & \int_{0}^{t}\int_{0}^{L}\left(-\mathbf{W}_{,t}^T\,\mathbf{U} - \mathbf{W}_{,z}^T\,\mathbf{F} + \mathbf{W}_{,z}^T\,\mathbf{K}\,\mathbf{U}_{,z}-\mathbf{W}^T\,\mathbf{G}\right)\,dz\,dt + \int_{0}^{T} \left[\mathbf{W}^T\left(\mathbf{F}-\mathbf{K}\mathbf{U}_{,z}\right)\right]_{0}^{L}\,dt + \\
   & \int_{0}^{L}\mathbf{W}^T(z,T)\mathbf{U}(z,T)\,dz - \\
   & \int_{0}^{L}\mathbf{W}^T(z,0)\,\mathbf{U}^0(z)\,dz = 0
   \end{eqnarray}


where the initial condition is given by :math:`\mathbf{U}^0(z) = \left[S^0(z),Q^0(z)\right]^T`. 
The boundary conditions are not specified at this point.

Disjoint Decomposition
----------------------

We adopt the disjoint decomposition approach described in 2.3 to derive appropriate outflow boundary conditions. First, we divide our spatial domain :math:`\Omega=[0,L]` into an upstream **numerical** domain :math:`\Omega^{n}: z\in(0,B)`, and a downstream **analytic** domain :math:`\Omega^{a}: z\in(B,L)`.

The boundary that separates these domains is defined as :math:`\Gamma_{B} : z = B`. We define a disjoint decomposition of our variables, for example for our unknown solution vector, :math:`\mathbf{U}`

.. math::

   \mathbf{U} = \mathbf{U}^{n} + \mathbf{U}^{a}

so that

.. math::

   \mathbf{U} = 
   \begin{cases}
   \mathbf{U}^{n} & z\in\Omega^{n}\\
   \mathbf{U}^{a} & z\in\Omega^{a}
   \end{cases}

We use a similar decomposition for our weighting function, :math:`\mathbf{W}` , and insert these expressions into our variational form. 

The disjoint nature of this expression is used to derive a new variational form for the 1D numerical domain: we obtain the original variational form specialized to the 1D numerical domain :math:`\Omega^{n}` with the addition of a boundary term accounting for the interface to the 1D analytic domain, :math:`\Omega^{a}`

.. math::
   :label: eq18

   \begin{eqnarray}
   & \int_{0}^{t}\int_{0}^{B}\left(-\mathbf{W}_{,t}^{n\,T}\,\mathbf{U}^{n} - \mathbf{W}_{,z}^{n\,T}\,\mathbf{F}(\mathbf{U}^{n}) + \mathbf{W}_{,z}^{n\,T}\,\mathbf{K}\,\mathbf{U}^{n}_{,z}-\mathbf{W}^{n\,T}\,\mathbf{G}(\mathbf{U}^{n})\right)\,dz\,dt \\
   & -\int_{0}^{B}\,\mathbf{W}^{n\,T}(z,T)\,\mathbf{U}^{n}(z,T)\,dz + 
   \int_{0}^{B}\mathbf{W}^{n\,T}(z,0)\,\mathbf{U}^{n}(z,0)\,dz + \\
   & \int_{0}^{T}\left[\mathbf{W}^{n\,T}\left(\mathbf{F}(\mathbf{U}^{n}) - \mathbf{K}\,\mathbf{U}^{n}_{,z}\right)\right]_{z=0}\,dt - 
   \int_{0}^{T}\left[\mathbf{W}^{a\,T}\left(\mathbf{F}(\mathbf{U}^{a}) - \mathbf{K}\,\mathbf{U}^{a}_{,z}\right)\right]_{z=B}\,dt = 0
   \end{eqnarray}


Now, we enforce the continuity of the weighting function at the interface:

.. math::
   
   \mathbf{W}^{a}\vert_{z=B} = \mathbf{W}^{n}\vert_{z=B}

and define the operators :math:`\mathbf{M}` and :math:`\mathbf{H}` on the :math:`\Omega^{a}` domain based on the model of the
downstream domain:

.. math::
   
   \left[\mathbf{M}(\mathbf{U}^{a})\right]_{z=B} = \left[\mathbf{F}(\mathbf{U}^{a}) - \mathbf{K}\mathbf{U}^{a}_{,z}\right]_{z=B}

:math:`\mathbf{M}` acts on the solution variables and :math:`\mathbf{H}` depends only on other terms like initial conditions, boundary conditions, and physical properties in the downstream domain. 

Finally, we enforce the continuity of the flux at the boundary:

.. math:: 
   :label: eq20

   \left[\mathbf{M}(\mathbf{U}^{n})\right]_{z=B} = \left[\mathbf{M}(\mathbf{U}^{a})\right]_{z=B} 
   
The final result is

.. math::

   \begin{eqnarray}
   & \int_{0}^{t}\int_{0}^{B}\left(-\mathbf{W}_{,t}^{n\,T}\,\mathbf{U}^{n} - \mathbf{W}_{,z}^{n\,T}\,\mathbf{F}(\mathbf{U}^{n}) + \mathbf{W}_{,z}^{n\,T}\,\mathbf{K}\,\mathbf{U}^{n}_{,z}-\mathbf{W}^{n\,T}\,\mathbf{G}(\mathbf{U}^{n})\right)\,dz\,dt \\
   & -\int_{0}^{B}\,\mathbf{W}^{n\,T}(z,T)\,\mathbf{U}^{n}(z,T)\,dz + 
   \int_{0}^{B}\mathbf{W}^{n\,T}(z,0)\,\mathbf{U}^{n}(z,0)\,dz + \\
   & \int_{0}^{T}\left[\mathbf{W}^{n\,T}\left(\mathbf{F}(\mathbf{U}^{n}) - \mathbf{K}\,\mathbf{U}^{n}_{,z}\right)\right]_{z=0}\,dt - 
   \int_{0}^{T}\left[\mathbf{W}^{n\,T}\left(\mathbf{M}(\mathbf{U}^{n}) + \mathbf{H}\right)\right]_{z=B}\,dt = 0
   \end{eqnarray}



We see that the solution in the numerical domain depends on the operators :math:`\mathbf{M}` and :math:`\mathbf{H}` defined by the mathematical model of the downstream domain but not the solution
variable, :math:`\mathbf{U}^{a}`, in the downstream domain.

The Map from the “DtN” Method
-----------------------------

The operators :math:`\mathbf{M}` and :math:`\mathbf{H}` are based on the mathematical model of the downstream domain using an approach based on the “Dirichlet-to-Neumann” method [5]_, [6]_, [7]_, [9]_.
The physics of the downstream domain depends on the upstream domain. Thus, an explicit solution on the downstream domain cannot be obtained. Instead, a relationship between the unknowns that incorporates all the information of the model, the *map*, is
derived. 
The DtN map is then inserted into the flux term previously described :eq:`eq20`, to derive the operators :math:`\mathbf{M}` and :math:`\mathbf{H}`. 
In practice, the contribution of the diffusive flux term :math:`\mathbf{K}\mathbf{U}_{,z}` is observed to be negligible in the boundary integral and is hence omitted in deriving an expression for :math:`\mathbf{M}` and :math:`\mathbf{H}` from equations :eq:`eq13` and :eq:`eq20`:

.. math::
   :label: eq23

   \begin{eqnarray}
   M_1(Q,S) + H_1 & = Q\\
   M_2(Q,S) + H_2 & = (1 + \delta)\,\frac{Q^2}{S} + \frac{1}{\rho}\int_{p_0}^{p} \tilde{S}(p,z,t)\,dp
   \end{eqnarray}

Note that the boundary conditions are not exact since, at a minimum, a linear approximation is employed in the downstream domain whereas a nonlinear model is used in the upstream domain.

Resistance (0D, constant in time)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can distinguish between instantaneous and memory cases. An example of an instantaneous map is when a simple proportional relationship is prescribed between pressure at time :math:`t` and flow at the same point in time that represents the resistance to flow
of the downstream domain, :math:`Q(B,t) = p(B,t)/R`. The resistance :math:`R` can be measured, taken from the literature or derived for Poiseuille flow (steady flow) or other models [82]. 
Then using equations :eq:`eq3`, :eq:`eq7`, :eq:`eq8`, and integrating the pressure term in :eq:`eq23`:

.. math::

   \begin{eqnarray}
   M_1(S) & = \frac{\tilde{p}(S,B)}{R},\quad H_1 = 0\\
   M_2(S) & = \frac{4}{3}\,\frac{M_1(S)^2}{S} + frac{4\,\sqrt{\pi}}{3}\,\frac{E\,h}{\rho}\,\sqrt{S},\quad H_2 = -\frac{4}{3\,\rho}\,E\,h\,\pi\,r^{0}(B)
   \end{eqnarray}

Windkessel RCR circuit model (0D, fully transient)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Flow and pressure are related by the following relationship

.. math::

   Q(B,t) = \left[Q(B,0) - \frac{p^0(B)}{R}\right]\exp(-\alpha\,t) + \frac{p(B,t)}{R} - \frac{1}{R^2\,C}\,\int_{0}^{t} p(B,\tau)\exp(-\alpha(t-\tau))\,d\tau

   \alpha = \frac{R + R_d}{R\,R_d\,C}

Then using equations :eq:`eq3`, :eq:`eq7` and :eq:`eq8`, and integrating the pressure term in :eq:`eq23`:

.. math::

   M_1(S) = \frac{\tilde{p}[S(B,t),B,t]}{R} - \frac{1}{R^2\,C}\,\int_{0}^{t}\tilde{p}[S(B,\tau),B,\tau]\,\exp(-\alpha(t-\tau))\,d\tau

   H_1 = \left[Q(B,0) - \frac{p^0(B)}{R}\right]\exp(-\alpha\,t)

   M_2(S) = \frac{4}{3}\,\frac{[M_1(S) + H_1]^2}{S} + \frac{4\,\sqrt{\pi}}{3}\,\frac{E\,h}{\rho}\,\sqrt{S}

   H_2 = -\frac{4}{3\,\rho}\,E\,h\,\pi\, r^{0}(B)


The flow rate at time :math:`t` depends on the entire history of the pressure represented by the time integral in the above equations.

Impedance (1D, periodic)
^^^^^^^^^^^^^^^^^^^^^^^^

Another example of a memory map is the impedance model: the downstream domain is approximated using linear wave propagation theory and we further assume that the solution is periodic in time. We can then derive

.. math::
   :label: eq27

   Q(B,t) = \frac{1}{T}\int_{t-\tau}^{t}\,p(B,\tau)\,y(B,t-\tau)\,d\tau


The flow rate at time :math:`t` depends on the history of the pressure over one period. Here :math:`y(B,t)` is the inverse Fourier transform of the admittance function 
The representation formula for the operators then reads, using equations :eq:`eq3`, :eq:`eq7`, :eq:`eq8`, :eq:`eq23`, :eq:`eq27`:

.. math::

   M_1(S) = \frac{1}{T}\,\int_{t-\tau}^{t}\,\tilde{p}\left[S(B,\tau),B\right]\,y(B,t-\tau)\,d\tau,\quad H_1=0

   M_2(S) = \frac{4}{3}\,\frac{M_1(S)^2}{S} + \frac{4\,\sqrt{\pi}}{3}\,\frac{E\,h}{\rho}\,\sqrt{S},\quad H_2 = -\frac{4}{3\,\rho}\,E\,h\,\pi\, r^{0}(B)

The flow rate at time :math:`t` depends on the history of the pressure over one cardiac cycle represented by the time integral in the above equations.

Wave in a tube (1D, fully transient)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another example of a memory map is the more general one-dimensional wave equation. The derivation of a minimally reflecting boundary condition for the one-dimensional non-linear equations using the wave equation for the downstream domain has been inspired by the work of Givoli, Grote and colleagues [7]_, [8]_, [9]_ on exact nonreflecting boundary conditions for the linear wave equation. 
For this latter case, we approximate the downstream domain using one-dimensional linear wave propagation theory but do not assume periodicity in time. 
As an example, in the case where the downstream domain is a single elastic vessel with length :math:`l` and wave speed :math:`c`, going from the boundary point :math:`B` to the far end point :math:`L`, we derived a map with the related Green’s function that relates cross-sectional area and its derivative at the inlet of a segment:

.. math::

  \frac{\partial S}{\partial z}(B,t) = -\frac{S(B,t)}{l} + \exp(\gamma\, t)\,\int_{0}^{t}\int_{B}^{L}\,\frac{\partial G}{\partial z}(B,t,z_0,t_0)\,f_B(z_0,t_0)\,dz_0\,dt_0 + \mathcal{H}(t)


Furthermore we integrate the balance of momentum equation in time to obtain:

.. math::

  Q(B,t) = -c^2\,\int_{0}^{t}\,\frac{\partial S}{\partial z}(B,t)\,\exp(2\gamma(t-t_0))\,dt_0 + Q^0(B)\exp(2\,\gamma\,t)

We can then derive a map between the flow rate and the cross-sectional area using (3.29) and (3.30):

.. math::

   Q(B,t) = 

   c^2\,\int_{0}^{t}\left[\frac{S(B,t^{*})}{l} - \exp(\gamma\,t^{*})\,\int_{0}^{t^{*}}\int_{B}^{L}\frac{\partial G}{\partial z}(B,t^{*},z_0,t_0)\,f_{B}(z_0,t_0)\,dz_0\,dt_0\right]\exp(2\,\gamma\,(t-t^{*}))\,dt^{*} + 

   c^2\,\int_{0}^{t}\mathcal{H}(t^{*})\exp(2\,\gamma\,(t-t^{*}))\,dt^{*} + Q^0(B)\exp(2\,\gamma\,t)

After integrating by parts in time, the derivatives that constitute :math:`f_B(z_0,t_0)`, and using the Green’s function :eq:`eq18`, the final map reads:

.. math::
   :label: eq32

   Q(B,t) = \frac{c^2}{l}\,\int_{0}^{t}\left[1 + \sum_{n=1}^{\infty}2\right]\,S(B,t^{*})\exp(2\,\gamma\,(t-t^{*}))\,dt^{*}

   -\left(\frac{c}{l}\right)^3\,\int_{0}^{t}\exp(\gamma\,(2\,t - t^{*}))\int_{0}^{t^{*}}\,S(B,t_0)\exp(-\gamma\,t_0)\left[\sum_{n=1}^{\infty}\frac{2\,n^2\,\pi^2}{\sqrt{\lambda_n}}\,\sin(c\,\sqrt{\lambda_n}(t^{*}-t_0))\right]\,dt_0\,dt^{*}

    + Q^0(B)\exp(2\,\gamma\,t) + \Theta\,\left[S^0(B), \dot{S}_0(B), S_L(t),\dot{S}_L(t),\ddot{S}_L(t)\right],

and

.. math::

   \Theta\,\left[S^0(B), \dot{S}_0(B), S_L(t),\dot{S}_L(t),\ddot{S}_L(t)\right] = - \left[\sum_{n=1}^{\infty}\frac{2\,c}{l\,\sqrt{\lambda_n}}\,\sin(c\,\sqrt{\lambda_n}\,t)\right]\exp(\gamma\,t)\,S(B,0) 

   + \left[\sum_{n=1}^{\infty}\,\frac{2\,l}{c\,n^2\,\pi^2\,\sqrt{\lambda_n}}\left(\gamma\,\sin(c\,\sqrt{\lambda}\,t)\right) + c\,\sqrt{\lambda_n}\,\left(\cos(c\,\sqrt{\lambda_n}\,t) - \exp(\gamma\,t)\right)\right]\exp(\gamma\,t)\,\dot{S}(B,0) + 

   -c^2\int_{0}^{t}\mathcal{H}(t^{*})\exp(2\,\gamma\,(t-t^{*}))\,dt^{*}


The operators for the wave boundary condition can now be derived using :eq:`eq23` and :eq:`eq32`, assuming as for the upstream numerical domain that the initial cross-sectional area is the same as the reference cross-sectional area:

.. math::

   Q(B,t) = M_1(S) + H_1,\, \gamma = \frac{N}{2\,S^{0}},\, \forall n \in \mathbb{N}_{>0},\,\lambda = \frac{n^2\,\pi^2}{l^2} - \frac{\gamma^2}{c^2}

   M_1(S) = \frac{c^2}{l}\int_{0}^{t}\left[1 + \sum_{n=1}^{\infty}\,2\right]\,S(B,t^{*})\exp\left[2\gamma(t-t^{*})\right]\,dt^{*} 

   - \left(\frac{c^2}{l}\right)^3\,\int_{0}^{t}\exp\left[\gamma(2t - t^{*})\right]\int_{0}^{t^{*}}S(B,t_0)\exp(-\gamma\,t_{0})\left[\sum_{n=1}^{\infty}\frac{2\,n^2\,\pi^2}{\sqrt{\lambda_n}}\,sin\left(t^{*} - t_{0}\right)\right]\,dt_0\,dt^{*}

   H_1 = Q^{0}(B)\exp\left(2\,\gamma\,t\right) + \Theta\left[S^0(B), \dot{S}^{0}(B), S_L(t), \dot{S}_L(t), \ddot{S}_L(t)\right]

   M_2(S) = \frac{4}{3}\frac{\left[M_1(S) + H_1\right]^2}{S} + \frac{4\,\sqrt{\pi}}{3}\frac{E\,h}{\rho}\sqrt{S}

   H_2 = -\frac{4}{3\,\rho}\,E\,h\,\pi\,r^0(B)



The flow rate is a function of pressure history and depends also on waves coming from the far end boundary conditions and the initial conditions everywhere in the downstream domain. 
For simplicity, we implemented the equation above assuming that the initial state corresponded to the static solution around which the wave equation is derived, with zero initial derivative of the cross-sectional area and a constant distant cross-sectional area.

The *DtN* map has now been derived for a variety of boundary conditions. The reader interested in the effect of a different boundary condition can follow the same approach to derive the corresponding map. 
In particular, this approach can be applied for complex lumped models of the coronary bed, and can also be performed very similarly for lumped-parameter heart models at the inlet of the numerical domain.

Finite Element Discretization
-----------------------------

We employ a stabilized space-time finite element method, known for its robustness, based on the Discontinuous Galerkin method in time. The procedure presented herein employs ideas developed in Hughes and Mallet [10]_ and Hughes, Franca and Hulbert [11]_. 
We previously [4]_ described a space-time method with flow rate, pressure and resistance boundary conditions that employed a  different strong form (non conservative).
Here we retained the same stabilization term. The present formulation accommodates more general inflow and outflow boundary conditions. We use shape functions that are piecewise constant in time and piecewise linear in space. 
Let :math:`\tilde{\mathcal{V}}` be the finite-dimensional approximation of :math:`\mathcal{V}` restricted to :math:`(0,B)\times(t_n,t_n+1)`. Thus, the weak form for slab :math:`n+1`, from :math:`t_n` to :math:`t_n+1` reads:

Find :math:`\mathbf{U^h}` in :math:`\mathbf{V^h}` such that :math:`\forall\mathbf{W^h}` in :math:`\mathbf{V^h}`.

.. math::

  \int_{t_n^{+}}^{t_{n+1}^{-}}\int_{0}^{B}\left[\mathbf{W}_{,t}^{T}\,\mathbf{U}^{\mathbf{h}} + 
  \mathbf{W}_{,z}^{T}\,\mathbf{F}(\mathbf{U}) - 
  \mathbf{W}_{,z}^{T}\,\mathbf{K}\,\mathbf{U}_{,z} + 
  \mathbf{W}^{T}\,\mathbf{G}\left(\mathbf{U}\right)
  \right]\,dz\,dt

  -\int_{0}^{B}\mathbf{W}^{T}\left(z,t_{n+1}^{-}\right)\,\mathbf{U}\left(z,t_{n+1}^{-}\right)\,dz
  + \int_{0}^{B}\,\mathbf{W}^{T}\left(z,t_{n+1}^{+}\right)\,\mathbf{U}\left(z,t_{n+1}^{-}\right)\,dz

  + \int_{t_{n}^{+}}^{t_{n+1}^{-}}\left\{\mathbf{W}\left[\mathbf{F}(\mathbf{U}) - \mathbf{K}\,\mathbf{U}_{,z}\right]\right\}_{z = 0}\,dt
  - \int_{t_{n}^{+}}^{t_{n+1}^{-}}\left\{\mathbf{W}\left[\mathbf{M}(\mathbf{U}) + \mathbf{H}\right]\right\}_{z = B}\,dt = 0


For simplicity, we have dropped the superscript :math:`h`. After discretization in time, (3.34) becomes (the superscript :math:`n+1` refers to time slab :math:`n+1`):

.. math::
   :label: eq35

   \Delta t_n\int_{0}^{B}\left[\mathbf{W}_{,z}^{T,n+1}\,\mathbf{F}^{n+1}(\mathbf{U^{n+1}}) - \mathbf{W}_{,z}^{T,n+1}\,\mathbf{K}\,\mathbf{U}_{,z}^{n+1} + \mathbf{W}^{T,n+1}\,\mathbf{G}^{n+1}(\mathbf{U}^{n+1})\right]\,dz

   - \int_{0}^{B}\,\mathbf{W}^{T,n+1}\left(\mathbf{U}^{n+1} - \mathbf{U}^{n}\right)\,dz
   + \Delta t_{n}\left[\mathbf{W}^{T,n+1}\left(\mathbf{F}^{n+1}(\mathbf{U}^{n+1}) - \mathbf{K}\,\mathbf{U}_{,z}^{n+1}\right)\right]_{z=0}

   - \int_{t_n^{+}}^{t_{n+1}^{-}}\left[\mathbf{W}^{T,n+1}\left(\mathbf{M}^{n+1}(\mathbf{U}) + \mathbf{H}^{n+1}\right)\right]_{z=B} = 0

The boundary term appears in the box in this last equation. 
Stabilization terms are added to :eq:`eq35`. To that end, we define the matrices :math:`\mathbf{A}` and :math:`\mathbf{C_A}`:

.. math::

   \mathbf{A} = 
   \begin{bmatrix}
   0 & 1\\
   -(1 + \delta)\,\left(\frac{U_2}{U_1}\right)^2 + \frac{U_1}{\rho}\frac{\partial\tilde{p}}{\partial S} & (1 + \delta)\,\frac{2\,U_2}{U_1}\\
   \end{bmatrix}   

   \mathbf{C_A} = 
   \begin{bmatrix}
   -\frac{\psi}{U_1} & 0\\
   f-\frac{1}{\rho}\,\frac{\partial\tilde{p}}{\partial z} & \frac{N}{U_1}\\
   \end{bmatrix}   

We also define the matrix differential operator:

.. math::

   \mathcal{L}(\mathbf{U}) = \mathbf{I}\,\frac{\partial}{\partial t} + \mathbf{A}(\mathbf{U})\,\frac{\partial}{\partial z} - \mathbf{K}\,\frac{\partial^2}{\partial z^2} - \mathbf{C_A}(\mathbf{U})


Note that :math:`\mathcal{L}(\mathbf{U})\,\mathbf{U}` is the residual of the advective form of the partial differential equation system. For the current case of a piecewise constant approximation in time and a piecewise linear approximation in space this simplifies to

.. math::

   \mathcal{L}(\mathbf{U})\,\mathbf{U} = \mathbf{A}(\mathbf{U})\,\mathbf{U}_{,z} - \mathbf{C_A}(\mathbf{U})\,\mathbf{U}

The stabilization term takes the form:

.. math::

   \Delta t_n\sum_{e}\int_{\Omega_e}\left(\mathcal{L}(\mathbf{U})^T\,\mathbf{W}\right)^T\boldsymbol{\tau}\,\mathcal{L}(\mathbf{U})\,\mathbf{U}\,dz

The summation ranges over the element interiors and :math:`\tau = \tau(\mathbf{U})` is the stabilization matrix defined by:

.. math::

   \boldsymbol{\tau} = \left[\frac{2}{\Delta t_n}\mathbf{I} + \frac{2}{h}\vert\mathbf{A}\vert + 3\,\left(\frac{2}{h}\right)^2\,\mathbf{K} + \vert\mathbf{C_A}\vert\right]^{-1}

Here, the absolute value of a 2x2 matrix B can be obtained from the Cayley-Hamilton theorem,

.. math::

   \vert\mathbf{B}\vert = \frac{\mathbf{B}^2 + \sqrt{det(\mathbf{B}^2)}\,\mathbf{I}}{\sqrt{tr(\mathbf{B}^2) + 2\,\sqrt{det(\mathbf{B}^2)}}}

Therefore the final variational problem is: find :math:`\mathbf{U}^{n+1}` such that :math:`\forall\,\mathbf{W}`:

.. math::

  \Delta t_n\int_{0}^{B}\left[\mathbf{W}_{,z}^{T,n+1}\,\mathbf{F}^{n+1}(\mathbf{U^{n+1}}) - \mathbf{W}_{,z}^{T,n+1}\,\mathbf{K}\,\mathbf{U}_{,z}^{n+1} + \mathbf{W}^{T,n+1}\,\mathbf{G}^{n+1}(\mathbf{U}^{n+1})\right]\,dz

  - \int_{0}^{B}\,\mathbf{W}^{T,n+1}\left(\mathbf{U}^{n+1} - \mathbf{U}^{n}\right)\,dz
  + \Delta t_{n}\left[\mathbf{W}^{T,n+1}\left(\mathbf{F}^{n+1}(\mathbf{U}^{n+1}) - \mathbf{K}\,\mathbf{U}_{,z}^{n+1}\right)\right]_{z=0}

  - \int_{t_n^{+}}^{t_{n+1}^{-}}\left[\mathbf{W}^{T,n+1}\left(\mathbf{M}^{n+1}(\mathbf{U}) + \mathbf{H}^{n+1}\right)\right]_{z=B} = 0

  + \Delta\,t_{n}\sum_{e}\int_{\Omega_e}\left(\mathbf{W}_{,z}^{T}\,\mathbf{A}^{n+1} - \mathbf{W}^{T}\,\mathbf{C}_{A}^{n+1}\right)\boldsymbol{\tau}\left(\mathbf{A}^{n+1}\mathbf{U}_{,z}^{n+1} - \mathbf{C}_{A}^{n+1}\,\mathbf{U}^{n+1}\right)\,dz = 0


Using piecewise linear shape functions in space :math:`N_A, A = 1,\dots, m` with :math:`m` the number of nodes, the global nodal residual is:

.. math::

  \mathbf{R}_{A} = \Delta\,t_{n}\,\int_{0}^{B}\,N_{A,z}\left(\mathbf{F}^{n+1}(\mathbf{U}^{n+1}) - \mathbf{K}\,\mathbf{U}_{,z}^{n+1}\right) + N_{A}\,\mathbf{G}^{n+1}(\mathbf{U^{n+1}})\,dz

  - \int_{0}^{B}\,N_{A}\left(\mathbf{U}^{n+1} - \mathbf{U}^{n}\right)\,dz
  + \Delta\,t_n\left[N_{A}\left(\mathbf{F}^{n+1}(\mathbf{U}^{n+1}) - \mathbf{K}\,\mathbf{U}_{,z}^{n+1}\right)\right]_{z=0}

  - \int_{t_{n}^{+}}^{t_{n+1}^{-}}\left[N_{A}\left(\mathbf{M}^{n+1}(\mathbf{U}) + \mathbf{H}^{n+1}\right)\right]_{z=B}\,dt

  + \Delta\,t_{n}\,\sum_{e}\,\int_{\Omega_e}\left(N_{A,z}\,\mathbf{A}^{n+1} - N_{A}\,\mathbf{C}_{A}^{n+1}\right)\,\boldsymbol{\tau}\,\left(\mathbf{A}^{n+1}\,\mathbf{U}_{,z}^{n+1} - \mathbf{C}_{A}^{n+1}\,\mathbf{U}^{n+1}\right)\,dz = 0

These nonlinear equations are then solved with a modified Newton-Raphson technique [4]_. At each iteration k+1 in the time step n+1, the non-linear loop consists of two steps:

Solve for the increment :math:`\Delta\mathbf{U}_{C}^{n+1,k+1}`:

.. math::

   \mathbf{K}_{AC}^{n+1,k}\,\Delta\mathbf{U}_{C}^{n+1,k+1} = -\mathbf{R}_{A}^{n+1,k},\quad\text{with}\quad\mathbf{K}_{AC}^{n+1,k} = \frac{\partial\mathbf{R}_{A}^{n+1,k}}{\partial\mathbf{U}_{C}},\quad A,C=1,\dots,m

Update the solution:

.. math::

   \mathbf{U}_{C}^{n+1,k+1} = \mathbf{U}_{C}^{n+1,k} + \Delta\mathbf{U}_{C}^{n+1,k+1}

The matrices :math:`\mathbf{A}`, :math:`\mathbf{C^A}`, :math:`\mathbf{C^F}` (recall :eq:`eq13`) and :math:`\tau` are frozen in the calculation of the tangent matrix:

.. math::

  \mathbf{K}_{AC} = \Delta\,t_n\,\int_{0}^{B}\,N_{A,z}\,\left(\frac{\partial\,\mathbf{F}^{n+1,k}}{\partial\,\mathbf{U}_{C}^{n+1,k}} - \mathbf{K}\,N_{C,z}\right) + N_A\,N_C\,\mathbf{C}_{F}^{n+1,k}\,dz

  - \int_{0}^{B}\,N_A\,N_C\,\mathbf{I}\,dz + 

  \Delta\,t_n\left[N_A\left(\frac{\partial\,\mathbf{F}^{n+1,k}}{\partial\,\mathbf{U}_{C}^{n+1,k}} - \mathbf{K}\,N_{C,z}\right)\right]_{z=0}

  - \int_{t_n^{+}}^{t_{n+1}^{-}}\,\left[N_A\,\frac{\partial\,\mathbf{M}^{n+1,k}}{\partial\,\mathbf{U}_{C}^{n+1,k}}\right]_{z=B}\,dt

  + \Delta\,t_{n}\,\sum_{e}\,\int_{\Omega_e}\left(N_{A,z}\,\mathbf{A}^{n+1} - N_{A}\,\mathbf{C}_{A}^{n+1}\right)\,\boldsymbol{\tau}\,\left(\mathbf{A}^{n+1}\,\mathbf{U}_{,z}^{n+1} - \mathbf{C}_{A}^{n+1}\,\mathbf{U}^{n+1}\right)\,dz = 0

After the residual converges to a chosen tolerance, the scheme is advanced in time to solve for a new time step, initialized with the solution at the previous time step.

**Remark 1** : In practice, the residuals and the tangent matrices are coded at the element level. The detailed finite element  residuals and tangent matrices are presented for reference in Appendix A for each boundary condition.

**Remark 2** : So far, we have presented the derivation for a single segment. At a connection of multiple segments, pressure continuity and conservation of mass are enforced using Lagrange multipliers. Pressure, cross-sectional area and flow rate boundary conditions are treated as essential boundary conditions. 
Both of these features are exactly the same as in Wan et al. [4]_ and are therefore not repeated here.

References
----------

.. [1] T.J.R. Hughes and J. Lubliner, **On the One-Dimensional Theory of Blood Flow in the Larger Vessels** , `Mathematical Biosciences` , 18(1-2) (1973), 161-170.
.. [2] T.J.R. Hughes, **A Study of the One-Dimensional Theory of Arterial Pulse Propagation**, 1974, `U.C. Berkeley`, Ph.D. Thesis.
.. [3] M.S. Olufsen, **Structured Tree Outflow Condition for Blood Flow in Larger Systemic Arteries** , `American Journal of Physiology` , 276 (1999), H257-268.
.. [4] J. Wan, B.N. Steele, S.A. Spicer, S. Strohband, G.R. Feijoo, T.J.R. Hughes and C.A. Taylor, **A One-Dimensional Finite Element Method for Simulation-Based Medical Planning for Cardiovascular Disease** , `Computer Methods in Biomechanics and Biomedical Engineering` , 5(3) (2002), 195-206.
.. [5] D. Givoli and J.B. Keller, **A Finite Element Method for Large Domains** , `Computer Methods in Applied Mechanics and Engineering` , 76(1) (1989), 41-66.
.. [6] J.B. Keller and D. Givoli, **Exact Non-Reflecting Boundary-Conditions** , `Journal of Computational Physics` , 82(1) (1989), 172-192.
.. [7] D. Givoli, **Numerical Methods for Problems in Infinite Domains**, 1992, Elsevier Science.
.. [8] M. Grote and J. Keller, **Exact Nonreflecting Boundary Conditions for the Time Dependent Wave Equation** , `SIAM Journal on Applied Mathematics` , 55(2) (1995), 280-297.
.. [9] I. Patlashenko, D. Givoli and P. Barbone, **Time-Stepping Schemes for Systems of Volterra Integro-Differential Equations** , `Computer Methods in Applied Mechanics and Engineering` , 190 (2001), 5691-5718.
.. [10] T.J.R. Hughes and M. Mallet, **A New Finite Element Formulation for Computational Fluid Dynamics: III. The Generalized Streamline Operator for Advective-Diffusive Systems** , `Computer Methods in Applied Mechanics and Engineering` , 58 (1986), 305-328.
.. [11] T.J.R. Hughes, L.P. Franca and G.M. Hulbert, **A New Finite Element Formulation for Computational Fluid Dynamics: VIII. The Galerkin/Least-Squares Method for Advective-Diffusive Equations** , `Computer Methods in Applied Mechanics and Engineering` , 73(2) (1989), 173-189.