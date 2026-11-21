# 6dim-lattice-covering
Supporting computational package for the paper "The Schuermann-Vallentin lattice L_6^c gives the thinnest lattice sphere covering in dimension six" by Alexey Garber

Date: November 20, 2026

Computational details

Section 7 of the paper  
[common/edges.py](common/edges.py) - Computes 17822 edges that do not violate the inertia bound for the standard simplex sigma0. See Lemma 7.1 of the paper.  
[02_step2_types/01_apex.py](02_step2_types/01_apex.py) - Computes 23 types of facets (or neighbor apexes) compatible with the 17822 edges. See Lemma 7.2 of the paper.  
[02_step2_types/02_types_ABCDE.py](02_step2_types/02_types_ABCDE.py) - Discards all but 5 types of facets using maxdet bound. See Lemma 7.3 of the paper.  
[02_step2_types/03_types_DE.py](02_step2_types/03_types_DE.py) - Discards Types D and E of facets DFS search. See Lemma 7.6 of the paper.

Section 8 of the paper  
[03_step3_typeA/claim1.py](03_step3_typeA/claim1.py) - Proves that if a unimodular simplex has a facet of Type A, then it has another facet of Type A uniquely defined by the neighbor. See Claim 8.1 of the paper.  
[03_step3_typeA/claim2.py](03_step3_typeA/claim2.py) - Proves that if a unimodular simplex has two facets of Type A from Claim 8.1, then it has another facet of Type A uniquely defined by the first two. See Claim 8.2 of the paper.  
[03_step3_typeA/claim3.py](03_step3_typeA/claim3.py) - Proves that if a unimodular simplex has four facets of Type A from Claims 8.1, 8.2, and Corollary 8.3, then the remaining three faces are of Types B or C. See Claim 8.4 of the paper.  
[03_step3_typeA/claim4.py](03_step3_typeA/claim4.py) - Proves that every simplex has a facet of Type A. See Claim 8.5 of the paper.

Section 9 of the paper  
[04_step4_cells/free_sum_cells.py](04_step4_cells/free_sum_cells.py) - Proves that one free sum of three triangles guaranteed by Type A facets uniquely defines the neighboring clusters forming other free sums. This allows us to use Sch\"urmann-Vallentin theorem about local optimality of $\mathsf{L}_6^c$. See Claim 9.1 of the paper.

Outputs from scripts from Sections 8 and 9 are in the corresponding folders. Scripts from Section 7 create json files used by further scripts.

Section 10 of the paper  
[exact_form/Q6c_exact.sage](exact_form/Q6c_exact.sage) - Contains SageMath verification of the exact expression for the quadratic form Q_6^c using computations in the field Q(x) where x is a root of a polynomial of degree 19 with integer coefficients. The script verifies that the given quadratic form has prescribed Delaunay triangulation and that it maximizes determinant on its secondary cone.
