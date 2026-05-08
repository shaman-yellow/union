

为确保推断结果的稳健性，我们制定了如下分析与质控流程：

1. 因果效应估算
    *   **模型选择**：首选**逆方差加权法（IVW）**作为因果效应估算的主分析方法。⟦mark$blue('根据 Cochran’s Q 检验结果动态选择模型：若存在显著异质性（P < 0.05），采用随机效应模型（Random effects model）；否则采用固定效应模型（Fixed effects model）')⟧。
    *   **多模型验证**：同步应用 MR-Egger 回归、加权中位数法（Weighted median）及基于模式的估算法（Mode-based estimates），通过多种统计假设的交叉验证评估效应的一致性。

2. 统计筛选与多维度质控
在完成计算后，我们对初始结果进行了严格筛选过滤：
    *   **效应显著性过滤**：基于 **IVW 模型**，以 ⟦mark$blue('{ (if(use_fdr) glue::glue("IVW {adjust_method} 校正后的 P 值") else "IVW P 值") } < {p_cutoff}')⟧ 判定为具有统计学意义的因果关联。
    *   **水平多效性控制 (Pleiotropy Test)**：为了排除遗传变异通过非暴露途径影响结局的风险，我们实施了双重监控：
        *   **MR-Egger 截距检验**：{ (if(require_no_pleio) glue::glue("要求其⟦mark$blue('截距项 P > {pleio_cutoff}，以确保不存在显著的常数项偏移')⟧。") else "评估其截距项以监控潜在的多效性偏倚。") }
        *   **MR-PRESSO 全局检验 (Global Test)**：利用 **`MRPRESSO`** ⟦pkgInfo('MRPRESSO')⟧ 进行异常值检测，⟦mark$blue('要求其 Global Test P > {pleio_cutoff}，从整体上排除由多效性异常值引起的估计偏差')⟧。
    *   **因果方向与异质性验证**：
        *   **方向性检验**：{ (if(require_steiger) "执行 ⟦mark$blue('Steiger 检验，仅保留显示正确因果方向的条目，排除逆向因果干扰')⟧。" else "通过 Steiger 检验辅助验证因果方向的一致性。") }
        *   **异质性过滤**：{ (if(require_heterogeneity) glue::glue("⟦mark$blue('仅保留 Cochran’s Q 检验 P > {heterogeneity_cutoff} 的结果')⟧，确保工具变量间的效应估计具有高度一致性。") else "利用 Cochran’s Q 统计量评估分析组合的异质性。") }


