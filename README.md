## A low-rank tensor completion method via Strassen-Ottaviani flattening

This package includes a MATLAB implementation of the algorithms (TCSO-NN and TCSO-WNN) presented in the paper "A low-rank tensor completion method via Strassen-Ottaviani flattening" by 
Tiantian He, Shenglong Hu and Zheng-Hai Huang.

If you have any questions about implementation, please contact Tiantian He (969795191@qq.com).

### Introduction:
The TCSO folder contains two algorithms:
- TCSO_NN.m: TCSO using the nuclear norm
- TCSO_WNN_ALM.m: TCSO using the weighted nuclear norm

We also provide two test demos:
- demo.m: tests random missing case 
- structural_demo.m: tests structural and mixed missing cases

Some notes have already been added to the code for easy reading.

For more detailed implementations, please refer to He et al. (2024) paper, which can be founded in the same directory.

### Reference:
If you find this code useful, please cite it in your publications as follows:

```
@article{he2024low,
  title={A Low-Rank Tensor Completion Method via Strassen--Ottaviani Flattening},
  author={He, Tiantian and Hu, Shenglong and Huang, Zheng-Hai},
  journal={SIAM J. Imaging Sci.},
  volume={17},
  number={4},
  pages={2242--2276},
  year={2024},
  publisher={SIAM}
}
```
