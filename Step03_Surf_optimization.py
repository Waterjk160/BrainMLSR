import numpy as np
from scipy.spatial import KDTree
from nibabel.freesurfer.io import read_geometry, write_geometry
from mne.transforms import apply_trans
import nibabel as nib
from scipy.ndimage import sobel
import argparse
import os



# def compute_laplacian(vertices, adjacency):
#     """计算拉普拉斯平滑项及其梯度"""
#     laplacian_energy = 0.0
#     gradient = np.zeros_like(vertices)
    
#     for i in range(len(vertices)):
#         neighbor_sum = np.sum(vertices[adjacency[i]], axis=0)
#         n_neighbors = len(adjacency[i])
        
#         # 均匀权重
#         delta_i = vertices[i] - neighbor_sum / n_neighbors
        
#         # 累加拉普拉斯能量
#         laplacian_energy += np.dot(delta_i, delta_i)
        
#         # 梯度
#         # gradient[i] = 2 * delta_i
#         gradient[i] = delta_i
    
#     return laplacian_energy, gradient

def compute_laplacian(vertices, adjacency):
    laplacian_energy = 0.0
    gradient = np.zeros_like(vertices)
    
    for i in range(len(vertices)):
        # 检查索引范围
        valid_indices = [idx for idx in adjacency[i] if 0 <= idx < len(vertices)]
        if valid_indices:
            neighbor_sum = np.sum(vertices[valid_indices], axis=0)
            n_neighbors = len(valid_indices)
            
            delta_i = vertices[i] - neighbor_sum / n_neighbors
            
            laplacian_energy += np.dot(delta_i, delta_i)
            gradient[i] = delta_i
        else:
            print(f"Vertex {i} has no valid neighbors.")
    
    return laplacian_energy, gradient

def compute_distance_error(vertices, target_coords, adjacency):
    """计算距离误差及其梯度"""
    distance_error = 0.0
    gradient = np.zeros_like(vertices)
    
    for i in range(len(vertices)):
        current_distance = np.linalg.norm(vertices[i] - target_coords[i])
        
        # 计算期望距离 d_0 (基于相邻顶点到目标表面的距离平均值)
        if len(adjacency[i]) > 0:
            d_0 = np.mean([np.linalg.norm(vertices[j] - target_coords[j]) for j in adjacency[i]])
        else:
            d_0 = current_distance
        
        direction = (vertices[i] - target_coords[i]) / current_distance if current_distance > 1e-6 else np.zeros(3)
        
        # 累加距离误差
        distance_error += (current_distance - d_0) ** 2 / 2
        
        # 梯度
        gradient[i] = (current_distance - d_0) * direction
    
    return distance_error, gradient

def compute_image_gradient_error(vertices, Torig, image_data, grad_x, grad_y, grad_z):
    """计算图像梯度"""
    image_gradient_energy = 0.0
    gradient = np.zeros_like(vertices)
    for i in range(len(vertices)):
        point_vox_coord = xyz_to_vox_coord_float(Torig, vertices[i])
        first_gradient_vox_value = trilinear_interpolation(image_data, point_vox_coord)
        image_gradient_energy += first_gradient_vox_value
        second_gradient_x = trilinear_interpolation(grad_x, point_vox_coord)
        second_gradient_y = trilinear_interpolation(grad_y, point_vox_coord)
        second_gradient_z = trilinear_interpolation(grad_z, point_vox_coord)

        # 梯度
        second_gradient_vox = np.array([second_gradient_x, second_gradient_y, second_gradient_z])
        gradient[i] = second_gradient_vox
        # gradient[i] = vox_to_xyz_coord(Torig, second_gradient_vox)
        # print(f"second_gradient_vox, gradient={second_gradient_vox}", flush=True)
        # print(f"gradient, gradient={gradient}", flush=True)
        # if i>10:
        #     break
    
    # print(f"compute_image_gradient_error, gradient={gradient}", flush=True)
    
    return image_gradient_energy, gradient

def trilinear_interpolation(brainmask_data, vox_coord):  
    x = vox_coord[0]
    y = vox_coord[1]
    z = vox_coord[2]
    # 确定整数坐标  
    x_floor, y_floor, z_floor = int(x), int(y), int(z)  
      
    # 确定周围的8个整数坐标点  
    coords = [  
        (x_floor,     y_floor,     z_floor),  
        (x_floor + 1, y_floor,     z_floor),  
        (x_floor,     y_floor + 1, z_floor),  
        (x_floor + 1, y_floor + 1, z_floor),  
        (x_floor,     y_floor,     z_floor + 1),  
        (x_floor + 1, y_floor,     z_floor + 1),  
        (x_floor,     y_floor + 1, z_floor + 1),  
        (x_floor + 1, y_floor + 1, z_floor + 1)  
    ]  
      
    # 获取这些点的体素值  
    values = [get_vox_value(brainmask_data, coord) for coord in coords] 

    # 插值计算  
    # 首先在z=z_floor平面上进行两次二维线性插值，得到两个中间值  
    # 然后在z=z_floor+1平面上进行两次二维线性插值，得到另外两个中间值  
    # 最后在这四个中间值上进行一维线性插值，得到目标点的体素值  
      
    # 在z=z_floor平面上的插值  
    val00 = values[0] + (x - x_floor) * (values[1] - values[0]) / 1  
    val01 = values[2] + (x - x_floor) * (values[3] - values[2]) / 1  
    val_z0 = val00 + (y - y_floor) * (val01 - val00) / 1  
      
    # 在z=z_floor+1平面上的插值  
    val10 = values[4] + (x - x_floor) * (values[5] - values[4]) / 1  
    val11 = values[6] + (x - x_floor) * (values[7] - values[6]) / 1  
    val_z1 = val10 + (y - y_floor) * (val11 - val10) / 1  
      
    # 在z方向上进行最终插值  
    result = val_z0 + (z - z_floor) * (val_z1 - val_z0) / 1  
      
    return result  

# 从体素坐标得到体素的值
def get_vox_value(brainmask_data, vox_coord):  
    return brainmask_data[vox_coord[0], vox_coord[1], vox_coord[2]]

def xyz_to_vox_coord_float(Torig, xyz):
    vox_coord = apply_trans(np.linalg.inv(Torig), xyz)
    return vox_coord

# 体素到XYZ
def vox_to_xyz_coord(Torig, vox):
    xyz = apply_trans(Torig, vox)
    return xyz


# def compute_gradient_vector_xyz(input_image_path):
#     # 加载MRI图像
#     mri_img = nib.load(input_image_path)
#     Torig = mri_img.header.get_vox2ras_tkr()  # 获取从体素坐标到物理坐标的转换矩阵
#     image_data = mri_img.get_fdata()
    
#     # 计算x, y, z方向上的梯度
#     grad_x = sobel(image_data, axis=0)  # x方向上的梯度
#     grad_y = sobel(image_data, axis=1)  # y方向上的梯度
#     grad_z = sobel(image_data, axis=2)  # z方向上的梯度
    
#     # 对每个点的方向向量进行归一化
#     grad_magnitude = np.sqrt(grad_x**2 + grad_y**2 + grad_z**2)
    
#     # 防止除以零，设置一个很小的值作为默认值
#     grad_magnitude[grad_magnitude == 0] = 1
    
#     grad_x /= grad_magnitude
#     grad_y /= grad_magnitude
#     grad_z /= grad_magnitude
    
#     return Torig, image_data, grad_x, grad_y, grad_z

def compute_gradient_vector_xyz(input_image_path):
    # 加载MRI图像
    mri_img = nib.load(input_image_path)
    
    # 使用 affine 替代 get_vox2ras_tkr
    Torig = mri_img.affine.copy()
    # Torig = mri_img.header.get_vox2ras_tkr()

    image_data = mri_img.get_fdata()
    
    # 计算x, y, z方向上的梯度
    grad_x = sobel(image_data, axis=0)  # x方向上的梯度
    grad_y = sobel(image_data, axis=1)  # y方向上的梯度
    grad_z = sobel(image_data, axis=2)  # z方向上的梯度
    
    # 对每个点的方向向量进行归一化
    grad_magnitude = np.sqrt(grad_x**2 + grad_y**2 + grad_z**2 + 1e-8)  # 避免除以零
    
    grad_x /= grad_magnitude
    grad_y /= grad_magnitude
    grad_z /= grad_magnitude
    
    return Torig, image_data, grad_x, grad_y, grad_z

# def bm_to_Torig_data(brainmask_file):
#     # brainmask到Torig和data
#     brainmask = nib.load(brainmask_file)
#     brainmask_data = brainmask.get_fdata()  
#     Torig = brainmask.header.get_vox2ras_tkr()
#     return Torig, brainmask_data

# def build_adjacency(num_vertices, faces):
#     """构建邻接矩阵"""
#     # num_vertices = np.max(faces) + 1
#     adjacency = [[] for _ in range(num_vertices)]
#     for face in faces:
#         # 确保face是一个长度为3的列表或元组
#         if len(face) != 3:
#             print(f"Invalid face: {face}", flush=True)
#             continue
#         for i in range(3):
#             if face[(i+1)%3] not in adjacency[face[i]]:
#                 adjacency[face[i]].append(face[(i+1)%3])
#             if face[(i+2)%3] not in adjacency[face[i]]:
#                 adjacency[face[i]].append(face[(i+2)%3])
#     return adjacency

def build_adjacency(num_vertices, faces):
    """构建邻接列表"""
    adjacency = [[] for _ in range(num_vertices)]
    
    for face in faces:
        # 确保face是一个长度为3的列表或元组
        if len(face) != 3:
            print(f"Invalid face: {face}")
            continue
        
        # 遍历face的每个顶点，添加邻居关系
        for i in range(3):
            v1 = face[i]
            v2 = face[(i + 1) % 3]
            v3 = face[(i + 2) % 3]
            
            # 检查索引范围
            if v1 < 0 or v1 >= num_vertices or v2 < 0 or v2 >= num_vertices or v3 < 0 or v3 >= num_vertices:
                print(f"Index out of range in face: {face}")
                continue
            
            # 添加边（无重复）
            if v2 not in adjacency[v1]:
                adjacency[v1].append(v2)
            if v3 not in adjacency[v1]:
                adjacency[v1].append(v3)
            
            if v1 not in adjacency[v2]:
                adjacency[v2].append(v1)
            if v3 not in adjacency[v2]:
                adjacency[v2].append(v3)
            
            if v1 not in adjacency[v3]:
                adjacency[v3].append(v1)
            if v2 not in adjacency[v3]:
                adjacency[v3].append(v2)
    
    return adjacency

def project_gradient_to_line_vectorized(grads, directions):
    """
    将一组梯度向量批量投影到对应的指定方向上。
    
    Parameters:
    - grads: 形状为 (N, 3) 的梯度向量数组，N 是顶点数量。
    - directions: 形状为 (N, 3) 的方向向量数组。
    
    Returns:
    - 投影后的梯度向量数组，形状为 (N, 3)。
    """
    # 计算单位方向向量
    norms = np.linalg.norm(directions, axis=1, keepdims=True)
    # 避免除以0，将接近0的值设置为1（这样不会改变原始梯度）
    norms[norms < 1e-6] = 1
    unit_directions = directions / norms
    
    # 计算投影长度
    projection_lengths = np.sum(grads * unit_directions, axis=1, keepdims=True)
    
    # 计算投影后的梯度
    projected_grads = projection_lengths * unit_directions
    return projected_grads


def read_abnormal_csv(csv_file):
    """从 CSV 文件读取 abnormal 标记"""
    df = pd.read_csv(csv_file)
    return df['is_abnormal'].values  # 返回异常标记的数组

# 原始的参数
# def gradient_descent(inner_coords, outer_coords, v_inner, v_outer, faces, input_image_path, \
#                     alpha_inner=3, alpha_outer=3,  beta_inner=1, beta_middle=1, beta_outer=1, gamma_inner=0.5, gamma_outer=0.5,\
#                     learning_rate=0.01, iterations=80, tol=1e-6):
def gradient_descent(inner_coords, outer_coords, v_inner, v_outer, faces, input_image_path, \
                    alpha_inner=3, alpha_outer=3,  beta_inner=1, beta_middle=1, beta_outer=1, gamma_inner=0.5, gamma_outer=0.5,\
                    learning_rate=0.01, iterations=30, tol=1e-6):
    """使用梯度下降优化顶点位置"""
    num_vertices = v_inner.shape[0]
    Torig, image_data, grad_x, grad_y, grad_z = compute_gradient_vector_xyz(input_image_path)
    # 构建邻接矩阵
    adjacency = build_adjacency(num_vertices, faces)
    
    print(f"alpha_inner={alpha_inner}, alpha_outer={alpha_outer},  beta_inner={beta_inner}, beta_middle={beta_middle}, beta_outer={beta_outer}\n, \
          gamma_inner={gamma_inner}, gamma_outer={gamma_outer}, learning_rate={learning_rate}, iterations={iterations}, tol={tol}")

    for iter in range(iterations):
        # 计算能量和梯度
        laplacian_energy_inner, laplacian_gradient_inner = compute_laplacian(v_inner, adjacency)
        laplacian_energy_outer, laplacian_gradient_outer = compute_laplacian(v_outer, adjacency)
        
        distance_error_outer, distance_gradient_outer = compute_distance_error(v_outer, outer_coords, adjacency)
        distance_error_middle, distance_gradient_middle_outer = compute_distance_error(v_outer, v_inner, adjacency)
        distance_error_middle, distance_gradient_middle_inner = compute_distance_error(v_inner, v_outer, adjacency)
        distance_error_inner, distance_gradient_inner = compute_distance_error(v_inner, inner_coords, adjacency)

        image_gradient_energy_outer, image_gradient_gradient_outer = compute_image_gradient_error(v_outer, Torig, image_data, grad_x, grad_y, grad_z)
        image_gradient_energy_inner, image_gradient_gradient_inner = compute_image_gradient_error(v_inner, Torig, image_data, grad_x, grad_y, grad_z)

        
        total_energy = alpha_inner * laplacian_energy_inner + alpha_outer * laplacian_energy_outer + \
                       beta_outer * distance_error_outer + beta_middle * distance_error_middle + beta_inner * distance_error_inner - \
                       gamma_outer * image_gradient_energy_outer - gamma_inner * image_gradient_energy_inner
        
        # 改成 Imgae gradinet +
        gradient_inner = alpha_inner * laplacian_gradient_inner + beta_inner * distance_gradient_inner + beta_middle * distance_gradient_middle_inner + gamma_inner * image_gradient_gradient_inner
        gradient_outer = alpha_outer * laplacian_gradient_outer + beta_outer * distance_gradient_outer + beta_middle * distance_gradient_middle_outer + gamma_outer * image_gradient_gradient_outer
        
        # # 打印每项能量
        # print(f"Iteration {iter}:")
        # print(f"  Laplacian Inner: {laplacian_energy_inner:.4f}")
        # print(f"  Laplacian Outer: {laplacian_energy_outer:.4f}")
        # print(f"  Distance Outer:  {distance_error_outer:.4f}")
        # print(f"  Distance Middle: {distance_error_middle:.4f}")
        # print(f"  Distance Inner:  {distance_error_inner:.4f}")
        # print(f"  Image Gradient Outer: {image_gradient_energy_outer:.4f}")
        # print(f"  Image Gradient Inner: {image_gradient_energy_inner:.4f}")
        # print(f"  Total Energy:    {total_energy:.4f}")
        # print("-" * 40)
        # print(f"Iteration {iter}: Total Energy = {total_energy}")
        # 对梯度进行投影处理
        deform_direction = outer_coords - inner_coords
        gradient_inner = project_gradient_to_line_vectorized(gradient_inner, deform_direction)
        gradient_outer = project_gradient_to_line_vectorized(gradient_outer, deform_direction)

        # 更新顶点位置
        new_v_inner = v_inner - learning_rate * gradient_inner
        new_v_outer = v_outer - learning_rate * gradient_outer
        
        # 通过check进行判断
        # 验证新位置的顺序性和共线性
        for i in range(len(new_v_inner)):
            if not check_collinearity_and_order_single(inner_coords[i], new_v_inner[i], new_v_outer[i], outer_coords[i]):
                # 如果不满足条件，恢复原来的坐标
                new_v_inner[i] = v_inner[i]
                new_v_outer[i] = v_outer[i]


        # 检查停止条件
        if np.linalg.norm(new_v_inner - v_inner) < tol and np.linalg.norm(new_v_outer - v_outer) < tol:
            print("Converged due to small change in vertices")
            break
            
        v_inner = new_v_inner
        v_outer = new_v_outer
    
    return v_inner, v_outer


def check_collinearity_and_order_single(inner_gm, inner_gr, outer_gr, outer_gm, tolerance=1e-3):
    """
    Check if a single set of points are collinear and correctly ordered.
    
    Parameters:
    - inner_gm, inner_gr, outer_gr, outer_gm: Coordinates of the point on each surface.
    - tolerance: Allowed deviation for considering vectors collinear.
    
    Returns:
    A boolean indicating whether the point meets the conditions.
    """
    # 计算向量
    vec1 = inner_gr - inner_gm
    vec2 = outer_gr - inner_gm
    vec3 = outer_gm - inner_gm
    
    # 检查共线性: 使用向量叉乘接近0来判断
    cross1 = np.linalg.norm(np.cross(vec1, vec2))
    cross2 = np.linalg.norm(np.cross(vec2, vec3))
    
    if cross1 > tolerance or cross2 > tolerance:
        return False

    # 检查方向一致性: 通过内积判断
    dot1 = np.dot(vec1, vec2)
    dot2 = np.dot(vec2, vec3)
    
    if not (dot2 >= dot1 and dot1 >= 0):
        return False
        
    return True

def main(white_surf, pial_surf, initial_hypointense_inner, initial_hypointense_outer, T2_gradient_image, output_file_final_inner, output_file_final_outer):
    # 加载数据
    inner_coords, faces = read_geometry(white_surf)
    outer_coords, _ = read_geometry(pial_surf)
    initial_v_inner, _,  = read_geometry(initial_hypointense_inner)
    initial_v_outer, _,  = read_geometry(initial_hypointense_outer)


    # 执行梯度下降优化
    optimized_v_inner, optimized_v_outer = gradient_descent(inner_coords, outer_coords, initial_v_inner, initial_v_outer, faces, T2_gradient_image)

    # 保存结果
    write_geometry(output_file_final_inner, optimized_v_inner, faces)
    write_geometry(output_file_final_outer, optimized_v_outer, faces)
    print(f"output_file_final_inner={output_file_final_inner}", flush=True)
    print(f"output_file_final_outer={output_file_final_outer}", flush=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process surface geometries and optimize using gradient descent.")
    parser.add_argument('white_surf', type=str, help='Path to the inner white matter surface file.')
    parser.add_argument('pial_surf', type=str, help='Path to the pial surface file.')
    parser.add_argument('initial_hypointense_inner', type=str, help='Path to the initial inner hypointense layer surf file.')
    parser.add_argument('initial_hypointense_outer', type=str, help='Path to the initial outer hypointense layer surf file.')
    parser.add_argument('T2_gradient_image', type=str, help='Path to the T2 gradient image file.')
    parser.add_argument('output_file_final_inner', type=str, help='Path to save the optimized inner final surf file.')
    parser.add_argument('output_file_final_outer', type=str, help='Path to save the optimized outer final surf file.')
    args = parser.parse_args()

    if not all(map(os.path.exists, [args.white_surf, args.pial_surf, args.initial_hypointense_inner, args.initial_hypointense_outer, args.T2_gradient_image])):
        missing_files = [f for f in ['white_surf', 'pial_surf', 'initial_hypointense_inner', 'initial_hypointense_outer', 'T2_gradient_image'] if not os.path.exists(getattr(args, f))]
        raise FileNotFoundError(f"The following files do not exist: {', '.join(missing_files)}")

    main(args.white_surf, args.pial_surf, args.initial_hypointense_inner, args.initial_hypointense_outer, args.T2_gradient_image, \
                                args.output_file_final_inner, args.output_file_final_outer)
