# wcurvint: 轻量开源曲线求交 C/C++ 库

wcurvint 是一个轻量、易用的开源曲线求交 C/C++ 库，专注于解决二维常用曲线的交点计算稳定性问题，提供简洁的 C 语言接口，适合个人项目、科研验证及中小型应用的快速集成。

## 📋 功能介绍和特点

### 核心功能

- 支持多种基础二维曲线的两两求交计算：
  - 直线 ↔ 直线
  - 直线 ↔ 圆弧
  - 直线 ↔ 贝塞尔曲线
  - 直线 ↔ 有理贝塞尔曲线
  - 圆弧 ↔ 圆弧
  - 圆弧 ↔ 贝塞尔曲线
  - 圆弧 ↔ 有理贝塞尔曲线
  - 贝塞尔曲线 ↔ 贝塞尔曲线
  - 贝塞尔曲线 ↔ 有理贝塞尔曲线
  - 有理贝塞尔曲线 ↔ 有理贝塞尔曲线
- 丰富的求交结果信息：
  - 结果数量
  - 曲线是否奇异（包含导数为0向量的点）
  - 交点和重合段
  - 交点是交叉点还是采样点（供上层应用选择不同策略）

### 核心特点
- 强大的健壮性
  - 基础区间迭代防止丢解
  - 区间代数多项式消元防止方程组化简导致精度含义变化
  - 端点处特殊处理，保证拓扑稳定
  - 求交结果优化，求交结果更符合预期
- 最基础的 C 语言接口，方便 C/C++、Python、C# 等其他语言调用

## 💡 专业版介绍

该开源版本能够满足一般的曲线求交需求。
如果有更高的效率要求或跨平台需求，我们还提供了专业版的曲线求交库。
- 效率方面：<b>专业版的效率是开源版的5倍</b>。
- 跨平台方面：开源版只支持Windows平台，<b>专业版提供更多的平台支持</b>。

<b>专业版更多信息，请咨询邮箱：1179422870@qq.com</b>

## 🚀 快速开始

### 环境准备

1. 目前该项目只支持Windows平台，如果需要其他平台支持，请参考后面专业版说明。
2. 下载本仓库源码
3. 运行 build.sh 编译
4. 包含头文件 wcurvint.h 到你的项目中
5. 链接生成的 wcurvint.lib，运行时将 wcurvint.dll 和 wsolver.dll 放置在可执行文件同级目录

### 代码示例

<b>直线与圆弧求交</b>
```cpp
void test_line_line() {
    // 创建直线
    void* curve0 = new_line2d();
    set_line2d_start_point(curve0, 0, 0);
    set_line2d_end_point(curve0, 100, 100);
    //创建圆弧
    void* curve1 = new_arc2d();
    set_arc2d_center(curve1, 50, 50);
    set_arc2d_radius(curve1, 10);
    set_arc2d_start_angle(curve1, 0);
    set_arc2d_delta_angle(curve1, 6.28);
    //求交
    void* intersections = new_curve_curve_intersections(curve0, curve1, 1E-6);
    //打印求交结果
    int count = get_curve_curve_intersection_count(intersections);
    std::cout << "Intersection count: " << count << std::endl;
    for (int i = 0; i < count; ++i) {
        bool is_overlap;
        double start_t0;
        double start_t1;
        bool is_start_sample;
        double end_t0;
        double end_t1;
        bool is_end_sample;
        get_curve_curve_intersection(intersections, i, is_overlap,
            start_t0, start_t1, is_start_sample, end_t0, end_t1, is_end_sample);
        std::cout << "Intersection" << i << ": ";
        if (is_overlap) {
            std::cout << "[" << start_t0 << ", " << end_t0 << "]" << "[" << start_t1 << ", " << end_t1 << "]" << std::endl;
        }
        else {
            std::cout << "[" << start_t0 << "]" << "[" << start_t1 << "]" << std::endl;
        }
    }
    //释放数据
    free(intersections);
    free_curve2d(curve0);
    free_curve2d(curve1);
}
```

<b>贝塞尔曲线与贝塞尔曲线求交</b>
```cpp
void test_bezier_bezier() {
    //创建贝塞尔曲线
    void* curve0 = new_bezier_curve2d(3);
    set_bezier_curve2d_control_point(curve0, 0, 0, 0);
    set_bezier_curve2d_control_point(curve0, 1, 100, 100);
    set_bezier_curve2d_control_point(curve0, 2, 200, 100);
    set_bezier_curve2d_control_point(curve0, 3, 300, 0);
    //创建贝塞尔曲线    
    void* curve1 = new_bezier_curve2d(3);
    set_bezier_curve2d_control_point(curve1, 0, 100, 0);
    set_bezier_curve2d_control_point(curve1, 1, 200, 100);
    set_bezier_curve2d_control_point(curve1, 2, 300, 100);
    set_bezier_curve2d_control_point(curve1, 3, 400, 0);
    //求交
    void* intersections = new_curve_curve_intersections(curve0, curve1, 1E-6);
    //打印求交结果
    int count = get_curve_curve_intersection_count(intersections);
    std::cout << "Intersection count: " << count << std::endl;
    for (int i = 0; i < count; ++i) {
        bool is_overlap;
        double start_t0;
        double start_t1;
        bool is_start_sample;
        double end_t0;
        double end_t1;
        bool is_end_sample;
        get_curve_curve_intersection(intersections, i, is_overlap,
            start_t0, start_t1, is_start_sample, end_t0, end_t1, is_end_sample);
        std::cout << "Intersection" << i << ": ";
        if (is_overlap) {
            std::cout << "[" << start_t0 << ", " << end_t0 << "]" << "[" << start_t1 << ", " << end_t1 << "]" << std::endl;
        }
        else {
            std::cout << "[" << start_t0 << "]" << "[" << start_t1 << "]" << std::endl;
        }
    }
    //释放数据
    free(intersections);
    free_curve2d(curve0);
    free_curve2d(curve1);
}
```

## 📞 联系我们

- 邮箱：1179422870@qq.com