

>For assignments implementations: [tbbbk/Scotty3D_Benky](https://github.com/tbbbk/Scotty3D_Benky)

[TOC]

# 1. Basic Math Review

**I didn't write everything in detail for this section.**

## **1.1 Linear Algebra Review**

### 1.1.1 Norm

Which measures total size, length, volume, intensity, etc.

Warning: L2 Norm does not encode geometric length unless vectors are encoded in an orthonormal basis.

### 1.1.2 Linear Map

Key Idea:  linear maps take lines to lines while keeps the origin fixed

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717152040852.png" alt="image-20250717152040852" style="zoom:50%;" />

It doesn't matter whether we add the vectors or apply the linear map first.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717152109816.png" alt="image-20250717152109816" style="zoom:50%;" />

### 1.1.3 Gram-Schmidt

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717154153328.png" alt="image-20250717154153328" style="zoom:50%;" />

Warning: for large number of vectors / nearly parallel vectors, this is not the best algorithm

## **1.2 Vector Calculus**

### 1.2.1 Matrix Representation of Cross Product

$$
\mathbf{u} := (u_1, u_2, u_3)\rightarrow 
\widehat{\mathbf{u}} := \begin{bmatrix} 0 & -u_3 & u_2 \\ u_3 & 0 & -u_1 \\ -u_2 & u_1 & 0 \end{bmatrix}\\
 \mathbf{u} \times \mathbf{v} = \widehat{\mathbf{u}} \mathbf{v} = \begin{bmatrix} 0 & -u_3 & u_2 \\ u_3 & 0 & -u_1 \\ -u_2 & u_1 & 0 \end{bmatrix} \begin{bmatrix} v_1 \\ v_2 \\ v_3 \end{bmatrix}
$$

### 1.2.2 Determinant of a Linear Map

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717161401840.png" alt="image-20250717161401840" style="zoom:50%;" />

### 1.2.3 Derivative as Best Linear Approximation

*Taylor series*:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717162139822.png" alt="image-20250717162139822" style="zoom:50%;" />

Replacing complicated functions with a linear (and sometimes quadratic) approximation is a powerful trick in graphics algorithms.

### 1.2.3 Gradients of Matrix-Valued Expressions

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717162939364.png" alt="image-20250717162939364" style="zoom:50%;" />

resource: [matrixcookbook.pdf](https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf)

### 1.2.4 Vector Fields

In general, a vector filed assigns a vector to each point in space, for example, we saw a gradient field:
$$
f(x,y)=x^2+y^2\\
\text{vector filed: }\nabla f(x,y)=(2x,2y)
$$

### 1.2.5 Divergence 散度

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717170702343.png" alt="image-20250717170702343" style="zoom:50%;" />

Commonly written as $\nabla · X$.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717164601687.png" alt="image-20250717164601687" style="zoom:50%;" />
$$
\nabla = \left( \frac{\partial}{\partial u_1}, \dots, \frac{\partial}{\partial u_n} \right)\\
X(\boldsymbol{u}) = \big( X_1(\boldsymbol{u}), \dots, X_n(\boldsymbol{u}) \big)\\
\nabla \cdot X := \sum_{i=1}^n \frac{\partial X_i}{\partial u_i}\\
\mathbb{R}^n\rightarrow \mathbb{R}
$$

### 1.2.6 Curl 旋度

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717171020355.png" alt="image-20250717171020355" style="zoom:50%;" />

Commonly written as $\nabla \times X$
$$
\nabla = \left( \frac{\partial}{\partial u_1}, \frac{\partial}{\partial u_2}, \frac{\partial}{\partial u_3} \right) \\
X(\boldsymbol{u}) = \big( X_1(\boldsymbol{u}), X_2(\boldsymbol{u}), X_3(\boldsymbol{u}) \big)\\
\nabla \times X := \begin{bmatrix} 
\dfrac{\partial X_3}{\partial u_2} - \dfrac{\partial X_2}{\partial u_3} \\[6pt]
\dfrac{\partial X_1}{\partial u_3} - \dfrac{\partial X_3}{\partial u_1} \\[6pt]
\dfrac{\partial X_2}{\partial u_1} - \dfrac{\partial X_1}{\partial u_2} 
\end{bmatrix}\\
\text{If we only consider two dimensions:}\\\
\nabla \times X := \frac{\partial X_2}{\partial u_1} - \frac{\partial X_1}{\partial u_2}
$$
<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717171032182.png" alt="image-20250717171032182" style="zoom:50%;" />

> Divergence of $X$ is the same as curl of 90-degree rotation of $X$
>
> <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717171235174.png" alt="image-20250717171235174" style="zoom:50%;" />

### 1.2.7 Hessian in Coordinates

Hessian is operator that gives us partial derivatives of the gradient.
$$
f(\mathbf{x}): \mathbb{R}^n \to \mathbb{R} \\
\nabla^2 f := \begin{bmatrix} 
\dfrac{\partial^2 f}{\partial x_1 \partial x_1} & \cdots & \dfrac{\partial^2 f}{\partial x_1 \partial x_n} \\[8pt]
\vdots & \ddots & \vdots \\[8pt]
\dfrac{\partial^2 f}{\partial x_n \partial x_1} & \cdots & \dfrac{\partial^2 f}{\partial x_n \partial x_n} 
\end{bmatrix}
$$

# 2. Rasterization

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717193824058.png" alt="image-20250717193824058" style="zoom:50%;" />

## **2.1 Rasterization 101: Drawing a Triangle**

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717193737840.png" alt="image-20250717193737840" style="zoom:50%;" />

Why triangle?

- Can approximate any shape
- always planar, well-defined normal
- easy to interpolate data at corners
- barycentric coordinates

Key reason: once everything is reduced to triangles, can focus on making an extremely well-optimized <u>pipeline</u> for drawing them

### 2.1.1 Pinhole camera

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717194036226.png" alt="image-20250717194036226" style="zoom:50%;" />

### 2.1.2 Computing triangle coverage

Key question: Which pixels does the triangle overlap?

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717194208562.png" alt="image-20250717194208562" style="zoom:50%;" />

### 2.1.3 Aliasing

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717194627973.png" alt="image-20250717194627973" style="zoom:50%;" />

The sampling process lose some information, leading the reconstruction result not exactly accurate. Or we say undersampling high-frequency signals results in <u>aliasing</u>.

### 2.1.4 SuperSampling

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717195519778.png" alt="image-20250717195519778" style="zoom:50%;" />

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717195532772.png" alt="image-20250717195532772" style="zoom:50%;" />

Split a pixel into $N \times N$ grids (e.g. $2\times2$, $4\times4$), or generate multiple random sampling points in each pixel. For each sampling points, execute the a complete pipeline and average these points to rasterize the pixel.

## **2.2 Spatial Transformations**

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717200120459.png" alt="image-20250717200120459" style="zoom:50%;" />

### 2.2.1 Rotation

Properties:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717200316632.png" alt="image-20250717200316632" style="zoom:50%;" />

Matrix Representations:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717200623424.png" alt="image-20250717200623424" style="zoom:50%;" />

For rotation matrix, the transpose matrix equals inverse matrix.
$$
R^TR=I
$$
<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717200924923.png" alt="image-20250717200924923" style="zoom:50%;" />

> Orthogonal Transformations
>
> We said the rotation inverse matrix equals transpose matrix, but not all $R^TR=I$ does not mean it is a rotation.
>
> When $R^TR=I$,
>
> - Rotations additionally **preserve** orientation: $det(R)>0$
> - <u>Reflections</u> reverse orientation: $det(R)<0$
>
> $$
> R = \begin{bmatrix} -1 & 0 \\ 0 & 1 \end{bmatrix}\\
> R^\top R = \begin{bmatrix} (-1)^2 & 0 \\ 0 & 1 \end{bmatrix} = I
> $$
>
> <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717201550035.png" alt="image-20250717201550035" style="zoom:50%;" />

### 2.2.2 Scaling

$$
f(\mathbf{u})=a\mathbf{u}
$$

Matrix Representations:
$$
D = \begin{bmatrix} a & 0 & 0 \\ 0 & a & 0 \\ 0 & 0 & a \end{bmatrix}, \quad \mathbf{u} = \begin{bmatrix} u_1 \\ u_2 \\ u_3 \end{bmatrix}\\
D \mathbf{u} = \begin{bmatrix} a & 0 & 0 \\ 0 & a & 0 \\ 0 & 0 & a \end{bmatrix} \begin{bmatrix} u_1 \\ u_2 \\ u_3 \end{bmatrix} = \begin{bmatrix} a u_1 \\ a u_2 \\ a u_3 \end{bmatrix} = a \mathbf{u}
$$

> Spectral Theorem:
>
> A symmetric matrix $A=A^T$ has
>
> - orthonormal eigenvectors $e_1,...,e_n\in\mathbb{R}^n$
> - real eigenvalues $\lambda_1,...,\lambda_n\in\mathbb{R}$
>
> $$
> Ae_i=\lambda_ie_i
> $$
>
> Hence, every symmetric matrix performs a non-uniform scaling  along some set of orthogonal axes.
>
> If $A$ is positive definite ($\lambda_i>0$), this scaling is positive

### 2.2.3 Shear

$$
f_{\mathbf{u},\mathbf{v}}=\mathbf{x}+\langle\mathbf{v},\mathbf{x}\rangle\mathbf{u}
$$

### 2.2.4 Composition

We can composite transformations:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717203617714.png" alt="image-20250717203617714" style="zoom:50%;" />

How do we decompose a linear transformation into pieces?

### 2.2.5 Polar & Singular Value Decomposition

Polar decomposition decomposes any matrix $A$ into orthogonal matrix $Q$ (rotation) and symmetric positive-semidefinite matrix $P$ (scaling).
$$
A = QP
$$
Since $P$ is symmetric, can take this further via the spectral decomposition  $P=VDV^T$ ($V$ orthogonal (eigenvectors),  $D$ diagonal (eigenvalue)):
$$
A&=QVDV^T\\
&=UDV^T\\
$$
Where $U$ and $V^T$ are rotation matrix and $D$ axis-aligned scaling matrix.

### 2.2.6 Interpolating Transformations—Polar

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717204616722.png" alt="image-20250717204616722" style="zoom:50%;" />

### 2.2.7 Translation!

$$
f_\mathbf{u}(\mathbf{x})=\mathbf{x}+\mathbf{u}
$$

**This transformation is NOT linear!**

- Additivity
  $$
  f_{\mathbf{u}}(\mathbf{x} + \mathbf{y}) = \mathbf{x} + \mathbf{y} + \mathbf{u}\\
  f_{\mathbf{u}}(\mathbf{x}) + f_{\mathbf{u}}(\mathbf{y}) = (\mathbf{x} + \mathbf{u}) + (\mathbf{y} + \mathbf{u}) = \mathbf{x} + \mathbf{y} + 2\mathbf{u}
  $$
  
- Homogeneity
  $$
  f_{\mathbf{u}}(a\mathbf{x}) = a\mathbf{x} + \mathbf{u}\\
  af_{\mathbf{u}}(\mathbf{x}) = a(\mathbf{x} + \mathbf{u}) = a\mathbf{x} + a\mathbf{u}
  $$
  

Translation is **AFFINE**, not linear. Hence we cannot composite it with other linear transformations.

So, how do we composite all the transformations?

### 2.2.8 Homogeneous Coordinates

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717205256190.png" alt="image-20250717205256190" style="zoom:50%;" /><img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717205316388.png" alt="image-20250717205316388" style="zoom:50%;" />

Hint: Every points along a ray represents the same point.

Consider a point $\mathbf{p} = (x, y)$, and the plane $z=1$ in 3D. Any three $\hat{\mathbf{p}}=(a,b,c)$ such that $(a/c,b/c)=(x,y)$ are <u>homogeneous coordinates</u> for $\mathbf{p}$

---

If we apply a translation to a 2D point $\mathbf{p}$, all the homogeneous coordinates $\hat{\mathbf{p}}$ looks like shear transformation.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717211258585.png" alt="image-20250717211258585" style="zoom:50%;" />

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717211313093.png" alt="image-20250717211313093" style="zoom:50%;" />

**But sheer is a linear transformation!**

Let a point $ \mathbf{p} = (p_1, p_2) $be translated by vector $\mathbf{u} = (u_1, u_2) $, resulting in:   $$ \mathbf{p}' = (p_1 + u_1, p_2 + u_2) $$      For homogeneous coordinates $\widehat{\mathbf{p}} = (c p_1, c p_2, c) $ ( $ c \neq 0 $ ), the translated coordinates become:   
$$
\widehat{\mathbf{p}}' = (c p_1 + c u_1, c p_2 + c u_2, c)
$$
Notice that we’re shifting $\widehat{\mathbf{p}}'$ by an amount  $c\mathbf{u}$ then  that’s  proportional to the distance  along the third axis—a shear.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717212046801.png" alt="image-20250717212046801" style="zoom:50%;" />

> Homogeneous coordinates can also used to distinguish the vectors and points
>
> For vector, the homogeneous coordinate is $0$, for point the homogeneous coordinate is $1$. Because we cannot translate a vector! So we need to set homogeneous coordinate to $0$ to ignore the translation.

### 2.2.9 Transformation Composition Order

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717212946732.png" alt="image-20250717212946732" style="zoom:50%;" />

### 2.2.10 Scene Graph

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250717212830584.png" alt="image-20250717212830584" style="zoom:50%;" />

## **2.3 3D Rotations and Complex Representations**

### 2.3.1 Degree of Freedom

We need three degrees of freedom to specify a rotation in 3D

Two determine the "axis" direction and one determine how much spinning around the "axis".

### 2.3.2 Commutativity of Rotation—3D

Order of rotation matters in 3D. Verify by yourself.

- Rotate 90° around Y, then 90° around Z, then 90° around X
- Rotate 90° around Z, then 90° around Y, then 90° around X 

### 2.3.3 Euler Angles—Gimbal Lock

We can use the Euler Angles to represent a 3D rotation like:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718143603218.png" alt="image-20250718143603218" style="zoom:50%;" />

But sometimes we will encounter the Gimbal Lock!

Define rotation matrices for $x, y, z$ axes:  
$$
R_x = \begin{bmatrix} 
1 & 0 & 0 \\
0 & \cos\theta_x & -\sin\theta_x \\
0 & \sin\theta_x & \cos\theta_x 
\end{bmatrix}, \quad 
R_y = \begin{bmatrix} 
\cos\theta_y & 0 & \sin\theta_y \\
0 & 1 & 0 \\
-\sin\theta_y & 0 & \cos\theta_y 
\end{bmatrix}, \quad 
R_z = \begin{bmatrix} 
\cos\theta_z & -\sin\theta_z & 0 \\
\sin\theta_z & \cos\theta_z & 0 \\
0 & 0 & 1 
\end{bmatrix}\\
R_x R_y R_z = \begin{bmatrix} 
\cos\theta_y \cos\theta_z & -\cos\theta_y \sin\theta_z & \sin\theta_y \\
\cos\theta_z \sin\theta_x \sin\theta_y + \cos\theta_x \sin\theta_z & \cos\theta_x \cos\theta_z - \sin\theta_x \sin\theta_y \sin\theta_z & -\cos\theta_y \sin\theta_x \\
-\cos\theta_x \cos\theta_z \sin\theta_y + \sin\theta_x \sin\theta_z & \cos\theta_z \sin\theta_x + \cos\theta_x \sin\theta_y \sin\theta_z & \cos\theta_x \cos\theta_y 
\end{bmatrix}
$$
For $\theta_y = \frac{\pi}{2}$ (so $\cos\theta_y = 0, \sin\theta_y = 1$), the matrix simplifies to:  
$$
R_x R_y R_z \Big|_{\theta_y = \pi/2} = \begin{bmatrix} 
0 & 0 & 1 \\
\cos\theta_z \sin\theta_x + \cos\theta_x \sin\theta_z & \cos\theta_x \cos\theta_z - \sin\theta_x \sin\theta_y \sin\theta_z & 0 \\
-\cos\theta_x \cos\theta_z + \sin\theta_x \sin\theta_z & \cos\theta_z \sin\theta_x + \cos\theta_x \sin\theta_y \sin\theta_z & 0 
\end{bmatrix}
$$
Now, no matter how we adjust $\theta_x$ and $\theta_z$, it can only rotate in on plane

### 2.3.4 Complex Analysis

First, change the way of your thinking!

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718144456464.png" alt="image-20250718144456464" style="zoom:50%;" />

Instead, imagine it's just a quarter-turn in CCW direction:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718144509367.png" alt="image-20250718144509367" style="zoom:50%;" />

And complex numbers are just 2-vectors with $1$ and $i$ bases:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718144642310.png" alt="image-20250718144642310" style="zoom:50%;" />

And the multiple operation is a little different, but the rest operations are the same:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718144727627.png" alt="image-20250718144727627" style="zoom:50%;" />
$$
z_1 := (r_1, \theta_1)\\
z_2 := (r_2, \theta_2)\\
z_1 z_2 = (r_1 r_2, \theta_1 + \theta_2)
$$
So:
$$
e^{\mathrm{i}\pi}+1=0\\
\text{Specialization of Euler's formula:}\\
e^{\mathrm{i}\theta}=cos(\theta)+\mathrm{i}sin(\theta)\\
\text{Can use to ``implement'' complex product:}\\
z_1=ae^{\mathrm{i}\theta},z_2=be^{\mathrm{i}\phi}\\
z_1z_2=abe^{\mathrm{i}(\theta+\phi)}
$$
*All the angle use <u>rad</u>.

### 2.3.5 Rotation: Matrices vs.Complex

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718145557037.png" alt="image-20250718145557037" style="zoom:50%;" />

### 2.3.6 Quaternions 四元数<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718150346893.png" alt="image-20250718150346893" style="zoom:25%;" />

Hamilton’s insight: in order to do 3D rotations in a way that  mimics complex numbers for 2D, actually need FOUR coordinates.

One real, and three imaginary:
$$
\mathbb{H}:=span(\{1,\mathrm{i},\mathrm{j},\mathrm{k}\})\\
q = a+b\mathrm{i}+c\mathrm{j}+d\mathrm{k}\in\mathbb{H}\\
\text{Quaternion product determined by:}\\
\mathrm{i}^2=\mathrm{j}^2=\mathrm{k}^2=\mathrm{i}\mathrm{j}\mathrm{k}=-1\\
i\,j = k,j\,k = i, k\,i = j,\\
j\,i = -k,k\,j = -i, i\,k = -j.
$$
*product no longer commutes!

To encode 3D points $x, y, z$ with quaternions, map them to: 
$$
(x, y, z) \mapsto 0 + xi + yj + zk
$$
 where $i, j, k$are imaginary units of quaternions.  

A quaternion can be represented as a pair:  
$$
(\underbrace{\text{scalar}}_{\mathbb{R}}, \underbrace{\text{vector}}_{\mathbb{R}^3}) \in \mathbb{H}
$$

$$
(a, \mathbf{u})(b, \mathbf{v}) = \big( ab - \mathbf{u} \cdot \mathbf{v},\, a\mathbf{v} + b\mathbf{u} + \mathbf{u} \times \mathbf{v} \big)\\
\mathbf{u}\mathbf{v} = \mathbf{u} \times \mathbf{v} - \mathbf{u} \cdot \mathbf{v}
$$

### 2.3.7 3D Rotations via Quaternions

Given axis $\mathbf{u}$, angle $θ$, quaternion $q$ representing rotation is

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718151448585.png" alt="image-20250718151448585" style="zoom:75%;" />

Here we consider the normal vector $\mathbf{x}$ and pure imaginary quaternions. And the rotation is $\mathbf{\bar{q}xq}$

## **2.4 Perspective Projection**

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718154207396.png" alt="image-20250718154207396" style="zoom:50%;" />

Distant objects appear smaller!

### 2.4.1 View Frustum

View frustum is region the camera can see:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718154513536.png" alt="image-20250718154513536" style="zoom:50%;" />

- Top / bottom / left / right planes correspond to four sides of the image
- Near / far planes correspond to closest/furthest thing we want to draw

### 2.4.2 Clipping

When objects are not visible to the camera / in view frustum, we clip it out!

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718154613976.png" alt="image-20250718154613976" style="zoom:50%;" />

*Also near/far clipping

### 2.4.3 Mapping Frustum to Unit Cube

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718155934753.png" alt="image-20250718155934753" style="zoom:50%;" />

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718155958131.png" alt="image-20250718155958131" style="zoom:50%;" />

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718160010583.png" alt="image-20250718160010583" style="zoom:50%;" />

Solve  $A\mathbf{x_i} = \mathbf{y_i}$ for unknown entries of $A$

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718160118025.png" alt="image-20250718160118025" style="zoom:75%;" />

While we take perspective projection into account, it is:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718160354823.png" alt="image-20250718160354823" style="zoom:80%;" />

> For derivation: [OpenGL Projection Matrix](https://www.songho.ca/opengl/gl_projectionmatrix.html)

> Warp Up:
>
> <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718162819109.png" alt="image-20250718162819109" style="zoom:50%;" />

## **2.5 Texture Mapping**

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718184524758.png" alt="image-20250718184524758" style="zoom:50%;" />

### 2.5.1 Barycentric Coordinates

Very useful for interpolation.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718163120828.png" alt="image-20250718163120828" style="zoom:50%;" />

You can also regard it as the area proportions.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718163204817.png" alt="image-20250718163204817" style="zoom:50%;" />

### 2.5.2 Perspective Correct Interpolation

Due to perspective projection (homogeneous divide), barycentric interpolation of values  on a triangle with different depths is not an affine function of screen XY coordinates.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718163741440.png" alt="image-20250718163741440" style="zoom:50%;" />

We want to interpolate attribute values linearly in <u>3D object space</u>, not image space.  If we compute barycentric coordinates using 2D (projected) coordinates,  leads to (derivative) discontinuity in interpolation where quad was split:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718163836457.png" alt="image-20250718163836457" style="zoom:67%;" />

How do we do?

1. Evaluate $Z:=1/z$ and $P:=\phi/z$ at each vertex, where $\phi$ is the attribute
2. Interpolate $Z$ and $P$ using standard (2D) barycentric coordinates
3. Divide interpolated $P$ by interpolated $Z$

> For a derivation, see [Microsoft Word - lowk_persp_interp_06.doc](https://www.comp.nus.edu.sg/~lowkl/publications/lowk_persp_interp_techrep.pdf)

### 2.5.3 Texture Coordinates

Texture coordinate define a mapping from  surface coordinates to points in texture domain. Often defined by linearly interpolating texture  coordinates at triangle vertices.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718184656978.png" alt="image-20250718184656978" style="zoom:50%;" />

Each vertex has a coordinate $(u,v)$ in texture space

| <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718184746083.png" alt="image-20250718184746083" style="zoom:50%;" /> | <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718184921694.png" alt="image-20250718184921694" style="zoom:50%;" /> |
| ------------------------------------------------------------ | ------------------------------------------------------------ |

### 2.5.4 Magnification vs. Minification

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718185119536.png" alt="image-20250718185119536" style="zoom:50%;" />

- Magnification: camera is very close to scene object
- Minification: scene objects are very far away

### 2.5.5 Bilinear Interpolation (Magnification)

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718185412767.png" alt="image-20250718185412767" style="zoom:50%;" />

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718185421074.png" alt="image-20250718185421074" style="zoom:50%;" />

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718185429273.png" alt="image-20250718185429273" style="zoom:50%;" />

### 2.5.6 MIP Map (Minification)

When a pixel on the screen covers many pixels of the texture, we can average texture values of these pixels, which is expensive to compute. So, we can precompute it and choose one demanded.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718190043933.png" alt="image-20250718190043933" style="zoom:50%;" />

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718190111593.png" alt="image-20250718190111593" style="zoom:50%;" />

To get the mipmap level, we need to compute differences between texture coordinate values at neighboring samples.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718190406826.png" alt="image-20250718190406826" style="zoom:50%;" />

### 2.5.7 Trilinear Interpolation for MIP Map

We can first bilinear interpolate with mipmap level, and then interpolate between two different levels.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718191205744.png" alt="image-20250718191205744" style="zoom:50%;" />

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718191213947.png" alt="image-20250718191213947" style="zoom:50%;" />
$$
w = d-\lfloor d \rfloor
$$

### 2.5.8 Anisotropic Filtering

At grazing angles, samples may be stretched out by (very) different amounts along $u$ and $v$.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718191718477.png" alt="image-20250718191718477" style="zoom:50%;" />

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718191725258.png" alt="image-20250718191725258" style="zoom:50%;" />

We can sample multiple times along the longer direction and less times along the short direction.

## **2.6 Depth and Transparency**

### 2.6.1 Depth-Buffer

For each sample, depth-buffer stores the depth of the <u>closest</u> primitive seen so far.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718192513487.png" alt="image-20250718192513487" style="zoom:50%;" />

> There are also color buffer, which can also be used for super-sampling, like (4 samples per pixel)
>
> <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718192656644.png" alt="image-20250718192656644" style="zoom:50%;" />

### 2.6.2 Compositing

We can represent opacity as $\alpha$

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718192842741.png" alt="image-20250718192842741" style="zoom:50%;" />

An image can have a $\alpha$ channel:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718192939432.png" alt="image-20250718192939432" style="zoom:50%;" />

### 2.6.3 Fringing

| No fringing                                                  | Fringing                                                     |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718193040854.png" alt="image-20250718193040854" style="zoom:50%;" /> | <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718193055158.png" alt="image-20250718193055158" style="zoom:50%;" /> |

What cause this?

### 2.6.4 (Non-) Premultiplied Alpha

If we  Composite image $B$ with opacity $\alpha_B$ over image $A$ with opacity $\alpha_A$

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718193656480.png" alt="image-20250718193656480" style="zoom:50%;" />

| Non-premultiplied alpha                                      | Premultiplied alpha                                          |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718193559404.png" alt="image-20250718193559404" style="zoom:50%;" /> | <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718193609401.png" alt="image-20250718193609401" style="zoom:50%;" /><br /><img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718194726595.png" alt="image-20250718194726595" style="zoom:50%;" /><br /><img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718193637653.png" alt="image-20250718193637653" style="zoom:50%;" /> |

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718194402971.png" alt="image-20250718194402971" style="zoom:50%;" />

If we do not use the premulitplied alpha, there is fringe. Because the upsampled color mix the background color with original color. So we have to pre-multiplied color and then upsample it.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718194948309.png" alt="image-20250718194948309" style="zoom:50%;" />

> Premultiplied alpha is better!

## **2.7 Rasterization Pipeline Summary**

1.  Transform triangle vertices into camera space

   <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718200212341.png" alt="image-20250718200212341" style="zoom:50%;" />

2. Apply perspective projection transform to transform triangle vertices  into normalized coordinate space

   <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718200248013.png" alt="image-20250718200248013" style="zoom:50%;" />

3. Clipping

   <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718200322968.png" alt="image-20250718200322968" style="zoom:50%;" />

4. Transform to screen coordinates.  Perform homogeneous divide, transform vertex xy positions from  normalized coordinates into screen coordinates (based on screen w,h)

   <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718200346392.png" alt="image-20250718200346392" style="zoom:50%;" />

5. Setup triangle (triangle preprocessing). Before rasterizing triangle, can compute a bunch  of data that will be used by all fragments

   <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718200432487.png" alt="image-20250718200432487" style="zoom:50%;" />

6. Sample coverage and compute triangle color at sample point

   <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718200521161.png" alt="image-20250718200521161" style="zoom:50%;" />

   <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718200540346.png" alt="image-20250718200540346" style="zoom:50%;" />

7. Perform depth test (if enabled)

   <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718200604510.png" alt="image-20250718200604510" style="zoom:50%;" />

8.  Update color buffer* (if depth test passed) (\* Possibly using OVER operation for transparency)

   <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718200653638.png" alt="image-20250718200653638" style="zoom:50%;" />

> OpenGL/Direct3D graphics pipeline:
>
> <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718200731356.png" alt="image-20250718200731356" style="zoom:50%;" />

# 3. Geometry

## **3.1 Encode Geometry**

### 3.1.1 Implicit Representations

Points aren’t known directly, but satisfy some relationship. E.g., unit sphere is all points such that $x^2+y^2+z^2=1$.
$$
f(x,y,z)=2-1.23
$$
<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718202442048.png" alt="image-20250718202442048" style="zoom:50%;" />

Now, find a point on the plane. And we can observe that implicit surfaces make sampling hard.
$$
f(x,y,z) = x^2 + y^2 + z^2 - 1.
$$
<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718202548213.png" alt="image-20250718202548213" style="zoom:50%;" />

Now check if a point is inside or outside the unit sphere. And we can observe implicit surfaces make  inside/outside tests task easy.

**Common Implicit Representations:**

- Algebraic Surfaces

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718203411978.png" alt="image-20250718203411978" style="zoom:50%;" />

- Constructive Solid Geometry (Boolean operations)

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718203447219.png" alt="image-20250718203447219" style="zoom:50%;" />

- Blobby Surfaces (Gradually blend surfaces together)

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718203538248.png" alt="image-20250718203538248" style="zoom:50%;" />

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718203600446.png" alt="image-20250718203600446" style="zoom:50%;" />

- Blending Distance Functions

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718203730717.png" alt="image-20250718203730717" style="zoom:50%;" />

- Level Set Methods (Surface is found where interpolated values equal zero )

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718203756044.png" alt="image-20250718203756044" style="zoom:50%;" />

- Mandelbrot Set 

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718203847745.png" alt="image-20250718203847745" style="zoom:50%;" />

Pros: 

- description can be very compact (e.g., a polynomial) 
- easy to determine if a point is in our shape (just plug it in!) 
- other queries may also be easy (e.g., distance to surface) 
- for simple shapes, exact description/no sampling error 
- easy to handle changes in topology (e.g., fluid)  

Cons: 

- expensive to find all points in the shape (e.g., for drawing) 
- very difficult to model complex shapes

### 3.1.2 Explicit Representations

All points are given directly. E.g., points on sphere are $(cos(u)sin(v),sin(u)sin(v),cos(v)),\text{ for }0\le u\lt 2\pi\text{ and  }0\le v \le\pi$

Many explicit representations in graphics. Like  triangle meshes,  polygon meshes, subdivision surfaces, NURBS  point clouds…

Unsimilar to implicit representations, explicit representations make sampling easy but inside/outside test hard.

**Common Explicit Representations:**

- Point Cloud

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718204044544.png" alt="image-20250718204044544" style="zoom:50%;" />

- Polygon Mesh ( Store vertices and polygons)

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250718204111958.png" alt="image-20250718204111958" style="zoom:50%;" />

- Bézier Curves/Surfaces

  Bernstein Basis:

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719152619267.png" alt="image-20250719152619267" style="zoom:50%;" />

  It can use to interpolate different points:

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719152734279.png" alt="image-20250719152734279" style="zoom:50%;" />

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719152750819.png" alt="image-20250719152750819" style="zoom:50%;" />

  We can piece together many Bézier curves to interpolate lots of points (because High-degree Bernstein polynomials don’t interpolate well):

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719153052146.png" alt="image-20250719153052146" style="zoom:50%;" />

  Bézier Patches is Bézier patch is sum of (tensor) products of Bernstein bases:

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719153332361.png" alt="image-20250719153332361" style="zoom:50%;" />

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719153400985.png" alt="image-20250719153400985" style="zoom:50%;" />

  By connecting Bézier curves, can connect Bézier patches  to get a surface:

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719153622651.png" alt="image-20250719153622651" style="zoom:50%;" />

  > Basically, it is weight-average points (2D, 3D)

- Rational B-Splines

  Bézier can’t exactly represent conics—not even the circle!

  Solution: interpolate in homogeneous coordinates, then  project back to the plane:

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719154824673.png" alt="image-20250719154824673" style="zoom:50%;" />

- NURBS:  (N)on-(U)niform (R)ational (B)-(S)pline

  - knots at arbitrary locations (non-uniform) 
  - expressed in homogeneous coordinates (rational) 
  - piecewise polynomial curve (B-Spline) 

  w is homogeneous coordinate, controlling "strength" of a vertex

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719154946255.png" alt="image-20250719154946255" style="zoom:50%;" />

  We can use tensor product to the NURBS curve

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719155033869.png" alt="image-20250719155033869" style="zoom:50%;" />

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719155041179.png" alt="image-20250719155041179" style="zoom:50%;" />

- Subdivision:

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719155241193.png" alt="image-20250719155241193" style="zoom:50%;" />


### 3.1.3 Summary

Some representations work better  than others—depends on the task!

## **3.2 Meshes and Manifolds**

### 3.2.1 Manifold Assumption

If you zoom in far enough, can draw a regular coordinate grid. (Very rough definition)

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719162325950.png" alt="image-20250719162325950" style="zoom:50%;" />

This is not manifold:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719162427742.png" alt="image-20250719162427742" style="zoom:50%;" />

Which of shapes are manifold?

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719162458536.png" alt="image-20250719162458536" style="zoom:50%;" />

Or, we can say: **A manifold polygon mesh has fans, not fins**

1. Every edge is contained in only two polygons (no “fins”)  
2. The polygons containing each vertex make a single “fan”

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719162549472.png" alt="image-20250719162549472" style="zoom:50%;" />

### 3.2.2 Why Do We Need Manifold?

1.  To make some assumptions about our geometry to keep data  structures/algorithms simple and efficient 
2. In many common cases, doesn’t fundamentally limit what  we can do with geometry

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719162726353.png" alt="image-20250719162726353" style="zoom:50%;" />

### 3.2.3 Halfedge Data Structure

*There are lots data structure can be used to store meshes, here we only talk about halfedge.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719162914367.png" alt="image-20250719162914367" style="zoom:80%;" />

Halfedge makes mesh traversal easy:

- Visit all vertices of a face

  ```c++
  Halfedge* h = f->halfedge;
  do {
      h = h->next;
      h->vertex...
  } while (h != f-> halfedge);
  ```

- Visit all neighbors of a vertex:

  ```c++
  Halfedge* h = v->halfedge;
  do {
      h = h->twin->next;
  } while (h != v-> halfedge);
  ```

- …

**Halfedge connectivity is always manifold:**

- Keep following `next`, and you’ll get faces.
- Keep following `twin` and you’ll get edges.
- Keep following `next->twin` and you’ll get vertices.

### 3.2.4 Halfedge Meshes Edition

- Edge Flip

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719163515986.png" alt="image-20250719163515986" style="zoom:50%;" />

- Edge Split

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719163539056.png" alt="image-20250719163539056" style="zoom:50%;" />

- Edge Collapse

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719163601355.png" alt="image-20250719163601355" style="zoom:50%;" />

- More…

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719163733181.png" alt="image-20250719163733181" style="zoom:70%;" />

## **3.3 Digital Geometry Processing**

### 3.3.1 Remeshing as Resampling

- Undersampling destroys features 
- Oversampling bad for performance

We need "good sampling"—"good" mesh. We need good approximation of original shape! Keep only elements that contribute information about shape. Add additional information where, e.g., curvature is large.

Vertices exactly on the surface doesn't mean it is a good approximation.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719164456367.png" alt="image-20250719164456367" style="zoom:50%;" />

(Some attributes, like normal is not good).

- One rule of thumb: triangle shape—**Delaunay**

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719164710338.png" alt="image-20250719164710338" style="zoom:50%;" />

  For any triangle in the decomposition, the interior of its circumcircle does not contain any other point.

- Another rule of thumb: regular vertex degree: Degree 6 for triangle mesh, 4 for quad mesh

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719165526606.png" alt="image-20250719165526606" style="zoom:50%;" />

### 3.3.2 Upsampling via (Catmull-Clark/Loop) Subdivision

There are lots of ways to do the subdivision.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719165641337.png" alt="image-2025071916564137" style="zoom:70%;" />

- **Catmull-Clark subdivision**

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719170216363.png" alt="image-20250719170216363" style="zoom:60%;" />

  > For more detailed tutorial: [Catmull-Clark Subdivision: The Basics – CodeItNow](https://www.rorydriscoll.com/2008/08/01/catmull-clark-subdivision-the-basics/)

- **Loop Subdivision**

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719170417585.png" alt="image-20250719170417585" style="zoom:50%;" />

  We can use edge operations to complete the subdivision:

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719170455161.png" alt="image-20250719170455161" style="zoom:50%;" />

  (Don’t forget to update vertex positions!)

### 3.3.3 Simplification via Edge Collapse

Basically, we assign each edge a cost, collapse the edge with least cost, and repeat until we reach the target.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719170542776.png" alt="image-20250719170542776" style="zoom:50%;" />

And we use **Quadric Error Metrics** to determine edge's cost.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719183316611.png" alt="image-20250719183316611" style="zoom:60%;" />

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719183516959.png" alt="image-20250719183516959" style="zoom:50%;" />

> For a derivation, see [Scotty3D_Benky/assignments/A2/simplify.md at main · tbbbk/Scotty3D_Benky](https://github.com/tbbbk/Scotty3D_Benky/blob/main/assignments/A2/simplify.md)

### 3.3.4 Isotropic Remeshing Algorithm

How to make triangles uniform shape & size?

Repeat four steps: 

- Split any edge over 4/3rds mean edge length
- Collapse any edge less than 4/5ths mean edge length
- Flip edges to improve vertex degree 
- Center vertices tangentially

## **3.5 Geometric Queries**

### 3.5.1 Ray Equation

$$
r(t)=\mathbf{o}+t\mathbf{d}
$$

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719190555553.png" alt="image-20250719190555553" style="zoom:50%;" />

### 3.5.2 Intersection

For **implicit surface intersection**, $f(r(t))=0$ and solve for $t$.

For **explicit surface intersection** (e.g. triangle), things become much harder and we do care about **performance**!

We will introduce Spatial Acceleration Data Structures!

## **3.6 Spatial Acceleration Data Structures**

What we care about most is the ray-triangle intersection!

### 3.6.1 Affine Map for Triangle

We can parameterize triangle given by vertices $\mathbf{p}_0,\mathbf{p}_1,\mathbf{p}_2$ using  barycentric coordinates:
$$
f(u,v)=(1-u-v)\mathbf{p}_0+u\mathbf{p}_1+v\mathbf{p}_2\\
f(u,v)=\mathbf{p}_0+u(\mathbf{p}_1-\mathbf{p}_0)+v(\mathbf{p}_2-\mathbf{p}_0)\\
$$
Now it's like: 

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719191531371.png" alt="image-20250719191531371" style="zoom:50%;" />

So now the ray-triangle intersection is like:
$$
\mathbf{p}_0+u(\mathbf{p}_1-\mathbf{p}_0)+v(\mathbf{p}_2-\mathbf{p}_0)=\mathbf{o}+t\mathbf{d}\\

   \begin{bmatrix}
     \mathbf{p}_1 - \mathbf{p}_0 &
     \mathbf{p}_2 - \mathbf{p}_0 &
     -\mathbf{d}
   \end{bmatrix}
   \begin{bmatrix}u\\v\\t\end{bmatrix}
   = \mathbf{o} - \mathbf{p}_0\\
   M \begin{bmatrix}u\\v\\t\end{bmatrix}
   = \mathbf{o} - \mathbf{p}_0
   \quad\Longrightarrow\quad
   \begin{bmatrix}u\\v\\t\end{bmatrix}
   = M^{-1}\,(\mathbf{o}-\mathbf{p}_0)\\
    u \ge 0,\quad v \ge 0,\quad u + v \le 1,\quad t \ge 0
$$
We only need to solve for $u,v,\text{ and }t$.

### 3.6.2 Bounding Box

We can pre-compute a bounding box around all primitives. If a ray intersect with a bounding box, then we test each primitives within this bounding box to avoid meaningless tradeoff.

Then use calculate the ray-axis-aligned box intersection:

The uniform equation is:
$$
\mathbf{N}^T(\mathbf{o}+t\mathbf{d})=c
$$
Solve for the $t$, e.g. intersection with $x_0$
$$
\mathbf{N}^T=[1\space0\space0]^T\text{ (we only care about x-axis)}\\
c=x_0\\
t=\frac{x_0-\mathbf{o}_{\mathbf{x}}}{\mathbf{d}_{\mathbf{x}}}
$$
More examples:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719193429378.png" alt="image-20250719193429378" style="zoom:50%;" />

### 3.6.3 Bounding Volume Hierarchy (BVH)

> BVH implementation assignment is really the pain in the ass…

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719194340940.png" alt="image-20250719194340940" style="zoom:70%;" />

How do we build the better BVH? We need a better partition.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719194912528.png" alt="image-20250719194912528" style="zoom:50%;" />

A good partitioning minimizes the <u>cost</u> of finding the closest  intersection of a ray with primitives in the node.
$$
C =C_{trav} +p_AN_AC_{isect} +p_BN_BC_{isect}
$$


- $C_{trav}$ is the cost of traversing an interior node (e.g., bounding box test)
- $C_A$ and $C_B$ are the costs of intersection with the resultant child subtrees 
-  $p_A$ and  $p_B$ are the probability a ray intersects the bbox of the child nodes A and B

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719195228703.png" alt="image-20250719195228703" style="zoom:50%;" />
$$
P(hitA|hitB)=\frac{S_A}{S_B}
$$
The pipeline about building BVH:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719195332186.png" alt="image-20250719195332186" style="zoom:66%;" />

Beside BVH, there are also lots data structure to accelerate:

- K-D tree

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719195642869.png" alt="image-20250719195642869" style="zoom:50%;" />

- Uniform grid

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719195654373.png" alt="image-20250719195654373" style="zoom:50%;" />

- Heuristic: Choose number of voxels ~ total number of primitives

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719195705443.png" alt="image-20250719195705443" style="zoom:50%;" />

- Quad-tree / octree

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719195733076.png" alt="image-20250719195733076" style="zoom:50%;" />

# 4. Ray Tracing

## **4.1 Color**

> For the color section, I literally omit a lot of contents… Because I didn't listen to the color lecture very carefully orz

Light is oscillating electric & magnetic field.

 KEY IDEA: frequency determines color of light

### 4.1.1 Intensity or Absorption

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719201526219.png" alt="image-20250719201526219" style="zoom:50%;" />

### 4.1.2 Emission and Reflection

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719201620893.png" alt="image-20250719201620893" style="zoom:50%;" />

### 4.1.3 Color Models

- RGB<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719201718777.png" alt="image-20250719201718777" style="zoom:50%;" />
- CMYK<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719201730038.png" alt="image-20250719201730038" style="zoom:50%;" />
- HSV<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719201749706.png" alt="image-20250719201749706" style="zoom:33%;" />
- SML
- XYZ
- …

### 4.1.4  Y’CbCr

- Y’ = luma: perceived luminance (same as L* in CIELAB)  
- Cb = blue-yellow deviation from gray
- Cr = red-cyan deviation from gray

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719201947336.png" alt="image-20250719201947336" style="zoom:50%;" />

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250719201936232.png" alt="image-20250719201936232" style="zoom:50%;" />

## **4.2 Radiometry**

### 4.2.1 Photon

Imagine every photon is a little **rubber ball** hitting the scene.

### 4.2.2 Radiant Energy $Q$

This it "the number of hits". Energy for **single** photon:
$$
Q=\frac{hc}{\lambda}\\
h\approx6.626\times 10^{-34}J·s\\
c\approx3.00\times 10^8m/s\\
\lambda\approx390-700\times 10^{-3}m\text{ (visible)}\\
\text{Unit: }\frac{(J\times s)(m/s)}{m}=J
$$
where: $h$ is Planck's constant, $c$ is speed of light, and $\lambda$ is wavelength (color!).

### 4.2.3 Radiant Flux $\Phi$ (Power) 

Energy per unit time (Watts) received by the sensor (or emitted by the light)
$$
\Phi = \lim_{\Delta t \to 0} \frac{\Delta Q}{\Delta t} = \frac{dQ}{dt}\\
Q = \int_{t_0}^{t_1} \Phi(t) \, \mathrm{d}t
$$

### 4.2.4 Irradiance $E$

Area density of radiant flux, given a senor of with area $A$, the average flux is :
$$
\frac{\Phi}{A}
$$
Irradiance ($E$) is given by taking limit of area at a single point on the sensor:
$$
E(p) = \lim_{\Delta \to 0} \frac{\Delta \Phi(p)}{\Delta A} = \frac{\mathrm{d}\Phi(p)}{\mathrm{d}A} \quad \left[ \frac{\mathrm{W}}{\mathrm{m}^2} \right]
$$

### 4.2.5 Lambert's Law

| <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720153803749.png" alt="image-20250720153803749" style="zoom:50%;" /> | <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720153821568.png" alt="image-20250720153821568" style="zoom:50%;" /> |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| $E=\frac{\Phi}{A}$                                           | $E=\frac{E}{A'}=\frac{\Phi cos\theta}{A}$                    |

### 4.2.6 Irradiance Falloff with Distance

Given flux $\Phi$:
$$
\begin{align}
E_1 &= \frac{\Phi}{4\pi r_1^2} \to \Phi = 4\pi r_1^2 E_1 \\
E_2 &= \frac{\Phi}{4\pi r_2^2} \to \Phi = 4\pi r_2^2 E_2 \\
\frac{E_2}{E_1} &= \frac{r_1^2}{r_2^2} = \left( \frac{r_1}{r_2} \right)^2
\end{align}
$$
Since same amount of energy is distributed  over larger and larger spheres, has to get darker  quadratically with distance.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720154233970.png" alt="image-20250720154233970" style="zoom:50%;" />

### 4.2.7 Solid Angles

| Radians                                                      | Steradians                                                   |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720154336701.png" alt="image-20250720154336701" style="zoom:50%;" /> | <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720154344424.png" alt="image-20250720154344424" style="zoom:50%;" /> |
| $\theta=\frac{l}{r}$                                         | $\Omega=\frac{A}{r^2}$                                       |

Differential solid angle:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720154646024.png" alt="image-20250720154646024" style="zoom:50%;" />
$$
\begin{align}
\mathrm{d}A&=(r\mathrm{d}\theta)(r\mathrm{sin}\theta\mathrm{d}\phi)\\
 & =r^2\mathrm{sin}\theta\mathrm{d}\theta\mathrm{d}\phi\\
 \mathrm{d}\omega&=\frac{\mathrm{d}A}{r^2}=\mathrm{sin}\theta\mathrm{d}\theta\mathrm{d}\phi\\
 
\Omega &= \int_{S^2} \mathrm{d}\omega \\
&= \int_{0}^{2\pi} \int_{0}^{\pi} \sin\theta \, \mathrm{d}\theta \, \mathrm{d}\phi \\
&= 4\pi
\end{align}
$$

### 4.2.8 Radiance $L$

Radiance is the solid angle density of irradiance:
$$
L(\mathrm{p}, \omega) = \lim_{\Delta \to 0} \frac{\Delta E_\omega(\mathrm{p})}{\Delta \omega} = \frac{\mathrm{d}E_\omega(\mathrm{p})}{\mathrm{d}\omega} \quad \left[ \frac{\mathrm{W}}{\mathrm{m}^2 \, \mathrm{sr}} \right]
$$
 where $E_\omega$ denotes that the differential surface area is  oriented to face in the direction $\omega$!

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720160237131.png" alt="image-20250720160237131" style="zoom:50%;" />

 In other words, radiance is energy along a ray defined by  origin point p and direction $\omega$!

Energy per unit time per unit area per unit solid angle…!

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720160824013.png" alt="image-20250720160824013" style="zoom:50%;" />

**Surface radiance**: 
$$
L(\mathrm{p}, \omega) = \frac{\mathrm{d}E(\mathrm{p})}{\mathrm{d}\omega \cos\theta} = \frac{\mathrm{d}^2 \Phi(\mathrm{p})}{\mathrm{d}A \, \mathrm{d}\omega \cos\theta}
$$
Reminder: Often need to distinguish between incident radiance and  exitant radiance functions at a point on a surface. In general: 
$$
L_i(\mathbf{p}, \omega)\ne L_o(\mathbf{p}\omega)
$$


<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720161616080.png" alt="image-20250720161616080" style="zoom:50%;" />



### 4.2.9 Spectral Radiance

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720161451051.png" alt="image-20250720161451051" style="zoom:50%;" />

Now, the radiance is radiant energy per unit time per unit area per unit solid angle. If we wanna get the **COLOR**, we need to add a "<u>per unit wavelength</u>"

### 4.2.10 Ambient Occlusion

 Assume spherical (vs. hemispherical) light source, “**at infinity**”. Irradiance is now **rotation, translation invariant**. Can pre-compute, “bake” into texture to enhance shading

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720161832092.png" alt="image-20250720161832092" style="zoom:50%;" />

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720161840990.png" alt="image-20250720161840990" style="zoom:50%;" />

### 4.2.11 Radiant Intensity $I$

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720162659268.png" alt="image-20250720162659268" style="zoom:50%;" />

Power per solid angle emanating from a point source.
$$
I(\omega) = \frac{\mathrm{d}\Phi}{\mathrm{d}\omega} \quad \left[ \frac{\mathrm{W}}{\mathrm{sr}} \right]
$$


## **4.3 The Rendering Equation**

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720165136981.png" alt="image-20250720165136981" style="zoom:70%;" />

### 4.3.1 Recursive Raytracing

Basic strategy: recursively evaluate rendering equation!

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720165211245.png" alt="image-20250720165211245" style="zoom:50%;" />

Renderer measures **radiance** along a ray:

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720165250920.png" alt="image-20250720165250920" style="zoom:50%;" />

### 4.3.2 Reflection

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720165404977.png" alt="image-20250720165404977" style="zoom:67%;" />

When the ray bounce in scene, how does the reflection of light affect the outgoing radiance?

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720165422869.png" alt="image-20250720165422869" style="zoom:50%;" />

What we are talking about is the scatter function in the rendering equation. Choice of reflection function determines surface appearance.

Some basic reflection functions:

| Reflection                                                   | Examples                                                     |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| Ideal specular: Perfect mirror                               | <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720165759154.png" alt="image-20250720165759154" style="zoom:50%;" /><img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720165806019.png" alt="image-20250720165806019" style="zoom:50%;" /> |
| Ideal diffuse: Uniform reflection in all directions          | <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720165811359.png" alt="image-20250720165811359" style="zoom:50%;" /><img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720165817245.png" alt="image-20250720165817245" style="zoom:50%;" /> |
| Glossy specular: Majority of light distributed in  reflection direction | <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720165824559.png" alt="image-20250720165824559" style="zoom:50%;" /><img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720165832597.png" alt="image-20250720165832597" style="zoom:50%;" /> |
| Retro-reflective: Reflects light back toward source          | <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720165838997.png" alt="image-20250720165838997" style="zoom:50%;" /><img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720165847571.png" alt="image-20250720165847571" style="zoom:50%;" /> |

### 4.3.3 Models of Scattering

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720170704604.png" alt="image-20250720170704604" style="zoom:80%;" />

What goes in must come out! (Total energy must be conserved) 

### 4.3.4 BRDF (Bidirectional Reflectance Distribution Function)

It encodes behavior of light that “bounces off” surface and calculate, when given incoming direction $ω_i$, how much light gets scattered in any given outgoing direction $ω_o$.
$$
\begin{gather}
f_r(\omega_i \to \omega_o) \geq 0 \\
\int_{\mathcal{H}^2} f_r(\omega_i \to \omega_o) \, \cos\theta \, d\omega_i \leq 1\\
\text{the sum ≤1 instead of =1 because the surface may absorb the energy }\\
\text{and convert it into heat or something.} \\
f_r(\omega_i \to \omega_o) = f_r(\omega_o \to \omega_i)
\end{gather}
$$
**Radiometric description of BRDF**:For a given change in the incident irradiance, how much does the exitant radiance change?

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720172507596.png" alt="image-20250720172507596" style="zoom:50%;" />
$$
f_r(\omega_i \to \omega_o) = \frac{\mathrm{d}L_o(\omega_o)}{\mathrm{d}E_i(\omega_i)} = \frac{\mathrm{d}L_o(\omega_o)}{\mathrm{d}L_i(\omega_i) \cos\theta_i} \quad \left[ \frac{1}{\mathrm{sr}} \right]
$$
**Common BRDF:**

- **Lambertian reflection** $f_r = \frac{\rho}{\pi}$

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720172946409.png" alt="image-20250720172946409" style="zoom:50%;" />

- **Specular reflection** $f_r(\theta_i, \phi_i; \theta_o, \phi_o) = \frac{\delta(\cos\theta_i - \cos\theta_o) \delta(\phi_i - \phi_o \pm \pi)}{\cos\theta_i}$

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720173012500.png" alt="image-20250720173012500" style="zoom:50%;" />

- **Refraction**:

  - Snell's Law $\eta_i \sin\theta_i = \eta_t \sin\theta_t$

    <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720173048433.png" alt="image-20250720173048433" style="zoom:50%;" />

  - Law of refraction: solve the $\mathbf{cos} \theta_t$ in the Snell's Law

  - Optical manhole: Only small “cone” visible, due to total internal reflection (TIR) (When light is moving from a more optically dense  medium to a less optically dense medium, light incident on boundary from large enough angle  will not exit medium.)

    <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720173310356.png" alt="image-20250720173310356" style="zoom:50%;" />

- **Fresnel reflection**: Many real materials:  reflectance increases w/  viewing angle

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720173356732.png" alt="image-20250720173356732" style="zoom:50%;" />

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720173341227.png" alt="image-20250720173341227" style="zoom:50%;" />

- Anisotropic reflection: Reflection depends on azimuthal angle $\phi$

  <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720173603656.png" alt="image-20250720173603656" style="zoom:50%;" />

### 4.3.5 Subsurface scattering

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720173654307.png" alt="image-20250720173654307" style="zoom:67%;" />

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720173706572.png" alt="image-20250720173706572" style="zoom:50%;" />

**BSSRDF**:
$$
L(x_o, \omega_o) = \int_{A} \int_{H^2} S(x_i, \omega_i, x_o, \omega_o) \, L_i(x_i, \omega_i) \, \cos\theta_i \, \mathrm{d}\omega_i \, \mathrm{d}A
$$

| ![image-20250720173822059](D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720173822059.png) | ![image-20250720173827694](D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720173827694.png) |
| ------------------------------------------------------------ | ------------------------------------------------------------ |

# 5. Optimization for Ray Tracing

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720174045210.png" alt="image-20250720174045210" style="zoom:50%;" />
$$
\color{red}{\text{How can we possibly evaluate this integral?}}
$$

## 5.1 Numerical Integration

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720174218971.png" alt="image-20250720174218971" style="zoom:50%;" />

Basic idea: 

- integral is “area under curve” 
- sample the function at many points 
- integral is approximated as  weighted sum

### 5.1.1 Gauss Quadrature

For any polynomial of degree n, we can always obtain the  exact integral by sampling at a special set of n points and  taking a special weighted combination.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720175349086.png" alt="image-20250720175349086" style="zoom:50%;" />

Weighted combination of sample points.

Key idea so far: To approximate an integral, we need 

1. **quadrature points**
2. **weights for each point**

$$
\int_{a}^{b} f(x) \, dx \approx \sum_{i=1}^{n} w_i f(x_i)
$$

### 5.1.2 Trapezoid rule

 Approximate integral of $f(x)$ by pretending function is piecewise affine.

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720203645901.png" alt="image-20250720203645901" style="zoom:67%;" />
$$
\begin{align}
h &= \frac{b - a}{n - 1} \\
\int_{a}^{b} f(x) \, dx &= h \left( \sum_{i=1}^{n-1} f(x_i) + \frac{1}{2} \left( f(x_0) + f(x_n) \right) \right)
\end{align}
$$
Work: $O(n)$

Error: $O(h^2)=O(\frac{1}{n^2})$

How about 2D?

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720203936134.png" alt="image-20250720203936134" style="zoom:60%;" />

Error is still O(h^2), but work now is $O(n^2)$ (n x n set of measurements).

How about k dimensions?

let $N=n^k$, then the error is $O(h^2)=O(\frac{1}{n^2})=O(\frac{1}{N^{\frac{2}{k}}})$

### 5.1.3 Curse of Dimensionality for Trapezoid Rule

How much does it cost to apply the trapezoid  rule as we go up in dimension? 

- 1D: $O(n) $
- 2D: $O(n^2)$
- … 
- kD: $O(n^k)$

For many problems in graphics (like  rendering), k is very, very big (e.g., tens or  hundreds or thousands). Applying trapezoid rule does not scale!

### 5.1.4 Sampling from Discrete Probability  Distributions

To randomly select an event, select $x_i$ if:
$$
P_{i-1} < \xi \leq P_i
$$
Here $\xi$ is a uniform random variable $\in [0,1)$

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720205359000.png" alt="image-20250720205359000" style="zoom:67%;" />

### 5.1.5 Sampling Continuous Random Variables  using the Inversion Method

Cumulative probability distribution function $P(x)=Pr(X<x)$, get the inverse function $P^{-1}(x)$.  $\xi$ is a uniform random variable $\in [0,1)$, solve for $x=P^{-1}(\xi)$.

We must know integral of $p(x)$ (to get $P(x)$), and also the inverse function

> First try: uniformly sampling unit circle
>
> <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720212617356.png" alt="image-20250720212617356" style="zoom:50%;" />
>
> <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720212632425.png" alt="image-20250720212632425" style="zoom:50%;" />
>
> <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720212645417.png" alt="image-20250720212645417" style="zoom:50%;" />
>
> *For the second line:![image-20250720212923119](D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720212923119.png)
>
> <img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720213015116.png" alt="image-20250720213015116" style="zoom:67%;" />

### 5.1.6 Rejection Sampling

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720213112564.png" alt="image-20250720213112564" style="zoom:67%;" />

<img src="D:\BingkuiTongPersonalWebsite\tbbbk.github.io\blogs\cg\images\image-20250720213143070.png" alt="image-20250720213143070" style="zoom:67%;" />

