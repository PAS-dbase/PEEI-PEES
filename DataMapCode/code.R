setwd("F:/课题/孙老师ECMO/新建文件夹")
#BiocManager::install("sf")
# 重新安装 gtable 包
# BiocManager::install("plotly")
library(readxl)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(plotly)
library(countrycode)
library(dplyr)
data=read_xlsx("area.xlsx")
colnames(data)
data$Normalized_Freq <- (data$Freq - min(data$Freq)) / (max(data$Freq) - min(data$Freq)) * 4 + 1
# 获取中国地图数据
world <- ne_countries(scale = "medium", returnclass = "sf")
china <- subset(world, name == "China")
china_provinces <- ne_states(country = "China", returnclass = "sf")
# 计算每个省份的医院数量
province_counts <- data %>%
  group_by(provinces) %>%
  summarise(total_patients = sum(Freq))
province_counts$name=c("Shanghai","Beijing","Sichuan","Anhui",
                       "Shandong","Guangdong","Guangxi","Jiangsu",
                       "Jiangxi","Hebei","Henan","Zhejiang","Hubei",
                       "Hunan","Gansu","Fujian","Liaoning","Chongqing",
                       "Shaanxi","Heilongjiang")
# 合并省份的病人数量与中国省份地图数据
china_provinces <- merge(china_provinces, province_counts, by = "name", all.x = TRUE)
china_provinces$total_patients


# 绘制中国地图，并根据省份病人数量变色，添加医院数据点
china_provinces$total_patients[which(is.na(china_provinces$total_patients))]=0
ggplot(data = china_provinces) +
  geom_sf(aes(fill = total_patients), color = "black") +  # 绘制中国地图，并根据病人数量变色
  scale_fill_gradient2(low = "white",high = "#CD5C5C", name = "病人数量") +  # 颜色渐变（病人少是蓝色，病人多是红色）
  # 绘制医院数据点
  geom_point(data = data,
             aes(x = 经度, y = 纬度, size = Freq, color = Freq), alpha = 0.7) + 
  scale_size_continuous(name = "病人数量", range = c(2, 10)) +  # 调整点的大小
  scale_color_gradient(low = "#92A8D1", high = "#CD5C5C") +  # 调整颜色（病人少是蓝色，病人多是红色）
  
  # 对病人数量大于40的医院标注名称
  geom_text(data = subset(data, Freq > 40), 
            aes(x = 经度, y = 纬度, label = 医院名称), 
            size = 3, vjust = -1, hjust = 0.5, color = "black") + 
  
  theme_minimal() +
  labs(title = "Data Distribution") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

# library(RColorBrewer) 
# # 绘制中国地图，并根据省份病人数量变色，添加医院柱状图
# ggplot(data = china_provinces) +
#   geom_sf(aes(fill = total_patients), color = "black") +  # 绘制中国地图，并根据病人数量变色
#   scale_fill_gradient2(low = "white",high = "#D1422F", name = "病人数量") +  # 颜色渐变（病人少是蓝色，病人多是红色）
#   # 绘制医院数据点
#   geom_segment(data = data, 
#                aes(x = 经度, xend = 经度+0.01, 
#                    y = 纬度, 
#                    yend = data$纬度+data$Freq*0.04,
#                    color ="#8BC8CB",size = Normalized_Freq)) +
#   # 对病人数量大于40的医院标注名称
#   geom_text(data = subset(data, Freq > 40), 
#             aes(x = 经度, y = 纬度, label = 医院编号), 
#             size = 3, vjust = -1, hjust = 0.5, color = "black") + 
#   
#   theme_minimal() +
#   labs(title = "Data Distribution") +
#   theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))



##################
#美化成柱状图   #
#################
# 自定义长三角形柱状图
create_triangle <- function(lon, lat, height, width = 0.5) {
  data.frame(
    x = c(lon - width / 2, lon + width / 2, lon),
    y = c(lat, lat, lat + height)
  )
}

triangles <- do.call(rbind, lapply(1:nrow(data), function(i) {
  triangle <- create_triangle(data$经度[i], data$纬度[i], data$Freq[i] / 40)
  cbind(triangle, name = data$医院编号[i], patients = data$Freq[i])}))

# # 计算柱状图顶部的位置，用于添加标签
# triangle_labels <- triangles %>%
#   group_by(name) %>%
#   summarize(
#     x = mean(x),  # 柱状图的中心位置
#     y = max(y),   # 柱状图的顶部位置
#     patients = first(patients)  # 患者数量
#   )

#triangles=unique(triangles)
library(ggnewscale)#允许使用多个比例尺

ggplot() +
  # 绘制中国地图
  geom_sf(data = china_provinces,aes(fill = total_patients), color = "black") +
  scale_fill_gradient(name = "Province Patients", low = "#F7BF95", high = "#D95B43", na.value = "white") +
  # 使用 ggnewscale 添加新的 fill 比例尺
  ggnewscale::new_scale_fill() +
  # 绘制长三角形柱状图
  geom_polygon(data = triangles, aes(x = x, y = y, group = name, fill = patients), alpha = 0.7) +
  # 设置颜色渐变
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  # 对病人数量大于40的医院标注名称
  geom_text(data = subset(data, Freq > 40),
            aes(x = 经度, y = 纬度, label = 医院编号),
            size = 3, vjust = -1, hjust = 0.5, color = "black") + 
  # 设置坐标范围和比例
  coord_sf(xlim = c(73, 135), ylim = c(18, 54)) +
  # 添加标题和标签
  labs(title = "Hospital and Province Patient Distribution in China",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  # 分别控制两个图例
  guides(fill = guide_legend(order = 1), fill = guide_legend(order = 2))
