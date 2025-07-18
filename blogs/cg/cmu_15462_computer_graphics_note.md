

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



