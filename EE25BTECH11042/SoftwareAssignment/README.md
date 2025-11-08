# Image Compression using Singular Value Decomposition (SVD)

This project implements **image compression** using **Singular Value Decomposition (SVD)** in **pure C**, without external libraries.  
The algorithm reads a **grayscale PGM image**, computes a **rank-reduced SVD approximation**, and reconstructs a compressed version of the image.

If we only keep the **top r singular values**, we get a **compressed approximation** of the image:


Algorithm Used: **Power Iteration SVD**

Instead of computing full SVD (which is computationally expensive), we approximate the top `r` singular vectors using **Power Iteration**:

This finds the **dominant eigenvector of (A^TA), which corresponds to the **right singular vector**.

Each new vector is **orthogonalized** to previously found vectors to obtain multiple singular components.

---

Program Workflow

1. Read a **PGM (P5)** grayscale image
2. Convert pixel data to a double precision matrix
3. For each rank component:
   - Start with a random vector (v)
   - Apply repeated **Power Iteration** to converge to a singular vector
   - Compute singular value
   - Compute matching left singular vector
4. Reconstruct compressed image:
5. Write the output back to **PGM** format

---

Advantages:
1. Pure C implementation
2. No external math/image libraries required
3. Works directly on numeric matrix form
4. Excellent control over compression level

Limitations
1. Works only for grayscale (PGM P5) images
2. Power iteration slows down for very large images (>1500×1500)
3. Not optimized with BLAS

References
1. Gilbert Strang, MIT Linear Algebra Lectures
2. SVD Theory in Image Compression
