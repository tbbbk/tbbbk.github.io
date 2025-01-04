import os
import shutil
import re

markdown_file = 'D:\\BingkuiTongPersonalWebsite\\tbbbk.github.io\\blogs\\cs229\\images\\CS229.md'
target_dir = r'D:\\BingkuiTongPersonalWebsite\\tbbbk.github.io\\blogs\\cs229\\images'

pattern = r'C:\\Users\\tbk\\AppData\\Roaming\\Typora\\typora-user-images\\([a-zA-Z0-9_-]+\.png)'

if not os.path.exists(target_dir):
    os.makedirs(target_dir)

with open(markdown_file, 'r', encoding='utf-8') as file:
    content = file.read()

image_files = re.findall(pattern, content)

for image_file in image_files:
    source_path = os.path.join(r'C:\\Users\\tbk\\AppData\\Roaming\\Typora\\typora-user-images', image_file)
    destination_path = os.path.join(target_dir, image_file)
    
    if os.path.exists(source_path):
        shutil.copy(source_path, destination_path)
        print(f'文件 {image_file} 已复制到 {destination_path}')
    else:
        print(f'源文件 {source_path} 不存在')

print('文件复制完成！')
