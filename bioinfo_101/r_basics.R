# 1. 最基本的赋值：在 R 里，我们用 <- 来赋值，而不是 =
# 创建一个向量，存5个基因的名字
gene_names <- c("Gene_A", "Gene_B", "Gene_C", "Gene_D", "Gene_E")

# 创建另一个向量，存这5个基因的表达量
expression_levels <- c(10.5, 200.1, 5.2, 88.9, 45.0)

# 2. 组装成“数据框”
my_data <- data.frame(
  Gene = gene_names,
  Expression = expression_levels
)

# 打印出来看看
print(my_data)

# 3. 生信最常见的操作：筛选数据
# 假设我们只想看表达量大于 20 的基因
high_expression_genes <- my_data[my_data$Expression > 20, ]

# 打印筛选后的结果
print(high_expression_genes)

