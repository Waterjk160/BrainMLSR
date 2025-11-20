# BrainMLSR
实现BrainMLSR流程
## 1. BrainMLSR流程

<p align="center">
  <img src="figures/pipeline.png" width="500" />
  <br>
  <em>Figure 1: pipeline.</em>
</p>

先得皮层内外表面。注意皮层内外表面的顶点必须一一对应并且具有相同的三角形面片关系。
然后再通过BrainMLSR进行多信号层重建

## 1. 皮层内外表面重建
为了便于广泛使用，我们以使用freesurfer为例，进行预处理并且得到皮层内外表面
现在有的图像，T1.nii.gz和T2FLAIR.nii.gz。
### 1.标准化到0.5mm分辨率

```
mri_convert T1.nii.gz T1_05.mgz -cs 0.5
mri_convert T2.nii.gz T2FLAIR_05.mgz -cs 0.5
```
### 2.T1配准到T2FLAIR图像
可以使用这里提供的配准函数。
```
python Register.py --fixed T2FLAIR_05.mgz --moving T1_05.nii.gz --output_dir T1_05_reg.nii.gz
```
### 3.用FreeSurfer对T1进行皮层重建
UII_5T替换成自己合适的，具体freesurfer重建皮层表面可以参考freesurfer官网。

```
recon-all -all -i T1_05_reg.nii.gz -s "UII_5T" -openmp 8
```

## 2.BrainMLSR

### 1. 初始表面提取
通过梯度方法先得到低信号层的初始内外表面。initial_surface_output_dir换成合适的输出路径。如果是lh的话，就对应把所有的rh改成lh
```
python Step01_Surf_Initialization.py \
    --white  rh.white  \
    --pial rh.pial \
    --T2flair T2FLAIR_05.mgz \
    --output_dir $initial_surface_output_dir \
    --hemisphere "rh"
```

### 2. 多信号层表面优化

先得到T2FLAIR的梯度图像
```
python $code_dir/Step02_GradImage.py T2FLAIR_05.mgz --input_path $T2_flair --output_path T2FLAIR_05_gradient.mgz
```

通过能量函数对曲面进行优化
```
python Step03_Surf_Initialization.py \
    --white_surf $INNER_WHITE --pial_surf $PIAL_SURF --initial_hypointense_inner $LAYER_45_WHITE --initial_hypointense_outer $LAYER_34_WHITE --T2_gradient_image $T2_GRADIENT_IMAGE
    --output_file_final_inner $lh_granular_inner --output_file_final_outer $lh_granular_outer
```
















