# REST2增强采样自动化项目框架

## 项目结构
```
REST2_AutoSampling/
├── config/
│   ├── default_config.yaml          # 默认参数配置
│   └── user_config.yaml             # 用户自定义配置
├── modules/
│   ├── __init__.py
│   ├── config_manager.py            # 配置管理模块
│   ├── structure_analyzer.py        # 结构分析模块
│   ├── solute_selector.py           # 溶质选择模块
│   ├── topology_builder.py          # 拓扑构建模块
│   ├── replica_generator.py         # 副本生成模块
│   ├── mdp_generator.py             # MDP文件生成模块
│   ├── gromacs_runner.py            # GROMACS执行模块
│   ├── temperature_controller.py    # 温度控制模块
│   └── analysis_tools.py            # 分析工具模块
├── templates/
│   ├── mdp_templates/               # MDP模板文件
│   │   ├── minimization.mdp
│   │   ├── equilibration.mdp
│   │   └── production.mdp
│   └── job_templates/               # 作业提交模板
├── scripts/
│   ├── main.py                      # 主执行脚本
│   ├── setup_rest2.py               # REST2设置脚本
│   └── monitor_simulation.py        # 模拟监控脚本
├── utils/
│   ├── __init__.py
│   ├── file_handler.py              # 文件处理工具
│   ├── logger.py                    # 日志系统
│   └── validator.py                 # 参数验证工具
└── examples/
    ├── peptide_example/
    └── small_molecule_example/
```

## 核心模块功能设计

### 1. 配置管理模块 (config_manager.py)
```python
class ConfigManager:
    def __init__(self, config_file=None):
        """初始化配置管理器"""
        pass
    
    def load_config(self, config_file):
        """加载配置文件"""
        pass
    
    def validate_config(self):
        """验证配置参数"""
        pass
    
    def get_parameter(self, key):
        """获取特定参数"""
        pass
```

**管理的关键参数：**
- `T_min`: 最低温度 (K)
- `T_max`: 最高温度 (K)
- `n_replicas`: 副本数目
- `replica_exchange_interval`: 交换间隔
- `cutoff_distance`: 选择氨基酸的cutoff距离
- `simulation_time`: 模拟时间
- `output_frequency`: 输出频率

### 2. 结构分析模块 (structure_analyzer.py)
```python
class StructureAnalyzer:
    def __init__(self, structure_file):
        """初始化结构分析器"""
        pass
    
    def identify_target_region(self, selection_criteria):
        """识别目标区域（peptide或小分子）"""
        pass
    
    def find_nearby_residues(self, target_atoms, cutoff):
        """找到cutoff范围内的氨基酸"""
        pass
    
    def generate_index_groups(self):
        """生成GROMACS索引组"""
        pass
```

### 3. 溶质选择模块 (solute_selector.py)
```python
class SoluteSelector:
    def __init__(self, structure_analyzer):
        """初始化溶质选择器"""
        pass
    
    def define_solute_region(self, target_selection, cutoff):
        """定义溶质区域"""
        pass
    
    def create_selection_groups(self):
        """创建选择组用于REST2"""
        pass
    
    def write_index_file(self, output_file):
        """写入索引文件"""
        pass
```

### 4. 副本生成模块 (replica_generator.py)
```python
class ReplicaGenerator:
    def __init__(self, config_manager):
        """初始化副本生成器"""
        pass
    
    def calculate_temperature_ladder(self):
        """计算温度梯度"""
        pass
    
    def setup_replica_directories(self):
        """设置副本目录结构"""
        pass
    
    def generate_topology_files(self):
        """为每个副本生成拓扑文件"""
        pass
```

### 5. 温度控制模块 (temperature_controller.py)
```python
class TemperatureController:
    def __init__(self, T_min, T_max, n_replicas):
        """初始化温度控制器"""
        pass
    
    def calculate_exponential_ladder(self):
        """计算指数分布温度梯度"""
        pass
    
    def calculate_scaling_factors(self):
        """计算REST2缩放因子"""
        pass
    
    def write_temperature_files(self):
        """写入温度相关文件"""
        pass
```

### 6. GROMACS执行模块 (gromacs_runner.py)
```python
class GromacsRunner:
    def __init__(self, config_manager):
        """初始化GROMACS执行器"""
        pass
    
    def run_preprocessing(self):
        """执行预处理步骤"""
        pass
    
    def run_minimization(self):
        """执行能量最小化"""
        pass
    
    def run_equilibration(self):
        """执行平衡模拟"""
        pass
    
    def run_rest2_production(self):
        """执行REST2生产模拟"""
        pass
```

## 配置文件示例 (default_config.yaml)

```yaml
# REST2 Configuration
rest2_settings:
  T_min: 300.0              # 最低温度 (K)
  T_max: 400.0              # 最高温度 (K)
  n_replicas: 8             # 副本数目
  replica_exchange_interval: 1000  # 交换间隔 (步数)
  
# 选择参数
selection:
  target_type: "peptide"    # "peptide" 或 "small_molecule"
  target_selection: "resname LIG"  # 目标分子选择语句
  cutoff_distance: 1.0      # cutoff距离 (nm)
  
# 模拟参数
simulation:
  minimization_steps: 5000
  equilibration_time: 1.0   # ns
  production_time: 100.0    # ns
  dt: 0.002                 # 时间步长 (ps)
  
# GROMACS设置
gromacs:
  gmx_mpi_command: "gmx_mpi"
  n_threads: 4
  gpu_usage: true
  
# 输出设置
output:
  trajectory_output_freq: 5000
  energy_output_freq: 1000
  log_output_freq: 1000
```

## 主执行流程 (main.py)

```python
#!/usr/bin/env python3

from modules.config_manager import ConfigManager
from modules.structure_analyzer import StructureAnalyzer
from modules.solute_selector import SoluteSelector
from modules.replica_generator import ReplicaGenerator
from modules.gromacs_runner import GromacsRunner

def main():
    # 1. 加载配置
    config = ConfigManager("config/user_config.yaml")
    
    # 2. 分析结构
    analyzer = StructureAnalyzer(config.get_parameter("structure_file"))
    
    # 3. 选择溶质
    selector = SoluteSelector(analyzer)
    selector.define_solute_region(
        config.get_parameter("target_selection"),
        config.get_parameter("cutoff_distance")
    )
    
    # 4. 生成副本
    replica_gen = ReplicaGenerator(config)
    replica_gen.setup_replica_directories()
    
    # 5. 执行模拟
    runner = GromacsRunner(config)
    runner.run_preprocessing()
    runner.run_minimization()
    runner.run_equilibration()
    runner.run_rest2_production()

if __name__ == "__main__":
    main()
```

## 使用示例

```bash
# 基本使用
python scripts/main.py --config config/user_config.yaml --structure input.gro --topology topol.top

# 设置REST2参数
python scripts/setup_rest2.py --T-min 300 --T-max 450 --replicas 12 --target "resname LIG" --cutoff 1.2

# 监控模拟进度
python scripts/monitor_simulation.py --run-dir ./rest2_simulation/
```

## 后续开发计划

1. **第一阶段**：实现核心配置管理和结构分析模块
2. **第二阶段**：开发溶质选择和温度控制功能
3. **第三阶段**：集成GROMACS执行和监控系统
4. **第四阶段**：添加分析工具和结果可视化
5. **第五阶段**：优化性能和错误处理

## 扩展性考虑

- 支持不同的力场选择
- 集成不同的增强采样方法
- 支持集群作业管理系统
- 添加实时监控和自动重启功能
- 支持多种输入格式 (PDB, GRO, etc.)
