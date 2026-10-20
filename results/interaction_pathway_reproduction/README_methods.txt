Models: primary ~ Donor + Celltype * Condition; sensitivity ~ Donor * Celltype + Condition * Celltype
msigdbr: 25.1.1
fgsea: 1.36.2
limma: 3.66.0
edgeR: 4.8.2
fgsea: minSize=15, maxSize=500, default eps; ranks sign(logFC)*sqrt(F), max signed rank per symbol.
CAMERA: voom, one feature per symbol by highest mean logCPM, estimated or fixed correlation 0.01.
IMPORTANT: inspect outputs and compare against manuscript; do not assume previously reported OXPHOS values match.
