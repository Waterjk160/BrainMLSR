import nibabel as nib
import numpy as np
import argparse
import os

def compute_gradient_magnitude(image_data):
    """
    计算给定3D图像数据的基础梯度大小（非Sobel）。
    
    :param image_data: 输入的3D图像数据
    :return: 梯度大小组成的3D数组
    """
    # 获取图像维度
    depth, height, width = image_data.shape
    
    # 初始化梯度数组，与输入图像大小相同
    grad_x = np.zeros_like(image_data, dtype=np.float32)
    grad_y = np.zeros_like(image_data, dtype=np.float32)
    grad_z = np.zeros_like(image_data, dtype=np.float32)
    
    # 计算x方向梯度 (depth轴) - 使用前向差分
    # 最后一层使用后向差分以保持边界
    grad_x[:-1, :, :] = image_data[1:, :, :] - image_data[:-1, :, :]
    grad_x[-1, :, :] = image_data[-1, :, :] - image_data[-2, :, :]
    
    # 计算y方向梯度 (height轴)
    grad_y[:, :-1, :] = image_data[:, 1:, :] - image_data[:, :-1, :]
    grad_y[:, -1, :] = image_data[:, -1, :] - image_data[:, -2, :]
    
    # 计算z方向梯度 (width轴)
    grad_z[:, :, :-1] = image_data[:, :, 1:] - image_data[:, :, :-1]
    grad_z[:, :, -1] = image_data[:, :, -1] - image_data[:, :, -2]
    
    # 计算梯度的平方和再开根号得到梯度大小
    gradient_magnitude = np.sqrt(grad_x**2 + grad_y**2 + grad_z**2)
    return gradient_magnitude

def save_image(output_path, affine, data):
    """
    保存图像数据为NIfTI或MGZ格式文件。
    
    :param output_path: 输出文件路径
    :param affine: 空间变换矩阵
    :param data: 要保存的图像数据
    """
    data = data.astype(np.float32)
    # 根据文件扩展名选择合适的图像类型
    if output_path.lower().endswith(('.nii', '.nii.gz')):
        img = nib.Nifti1Image(data, affine)
    else:
        img = nib.MGHImage(data, affine)
    nib.save(img, output_path)

def main(input_path, output_path):
    # 加载MRI图像
    mri_img = nib.load(input_path)
    mri_data = mri_img.get_fdata()
    
    # 计算基础梯度大小（非Sobel）
    gradient_magnitude = compute_gradient_magnitude(mri_data)
    
    # 保存结果
    save_image(output_path, mri_img.affine, gradient_magnitude)
    print(f"基础梯度图像已保存至 {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="处理MGZ或NIfTI文件，使用基础有限差分计算梯度大小。")
    parser.add_argument('--input_path', required=True, type=str, help='输入MGZ或NIfTI文件的路径。')
    parser.add_argument('--output_path', required=True, type=str, help='保存输出MGZ或NIfTI文件的路径。')
    args = parser.parse_args()
    
    if not os.path.exists(args.input_path):
        raise FileNotFoundError(f"输入文件不存在: {args.input_path}")
    
    main(args.input_path, args.output_path)
