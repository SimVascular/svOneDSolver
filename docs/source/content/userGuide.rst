User Guide and FAQs
================

What are nodes, segments and joints?
^^^^^^^^^^^^
* A *node* specifies the coordinates in (x,y,z) space of the connection between vessel segments. This corresponds to the location of the vessel centerline of the 3D model.
* A *segment* is a straight length of vasculature with tapered circular cross section. It connects different cross sections of the model and is symmetric across the xy plane.
* A *joint* is a connection between segments that ensures a unique pressure value and mass continuity. Joints can be within a single branch or at bifurcations of multiple branches. 


How to select a time step?
^^^^^
The size of your time step is influenced by several factors. For a straight, smooth, well-behaved model, a larger time step can be used. Additionally, for more stiff vessels (higher elastic modulus), a larger time step can be used. The more compliance is present in your model or the presence of a larger number of branches or very small branches will likely necessitate smaller time steps. Starting with a time step of 0.01 and decreaing by a factor of 1/2 until your model attains convergence is 


How to select the number element quadrature points?
^^^^^^^^^^
This number generally does not need to be changed and can remain at its default value of 4. 


How to select the number of elements per segment?
^^^^^^^^^^^^
For the most complete studies, a mesh convergence can be performed to ensure the proper number of elements per segment are used. If this is not desired, generally 10-25 elements per segment is sufficient. If you have very large segements, you will likely want to increase this to 50-150 elements. You can start at the smaller number of elements per segment and increase if you experience solver errors (see below).


Solver errors 
^^^^^^^^^^^^^^
If you receive an error while running a simulation, most commonly this can be solved by decreasing the time step (work in factors of 1/2) and increasing the number of elements per segment (work in factors of 2).

A common error, which can occur when there is a large difference in the inlet and outlet areas of a segment, is ``outlet areas going negative``. 

If this isn't working, then the geometry of the 1D model may need to be altered. This could involve adding additional segments to make the change in inlet and outlet area of the segments more gradual. It could also included truncating the model to remove branches or sections of branches with very small radii.
