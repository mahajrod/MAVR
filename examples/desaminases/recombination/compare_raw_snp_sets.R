# Test to check if previous(RUN2,3,4,5,6) and new data(RUN7) differ significantly(where it is possible)

PmCDA1_D6_old = c(1094, 785, 515, 1359)
PmCDA1_D6_RUN7 = c(292, 769, 205, 667, 274)

# Manna-Whitney test
wilcox.test(PmCDA1_D6_old, PmCDA1_D6_RUN7)
# W = 18, p-value = 0.06349

A1_D3_old = c(645, 606, 656, 642)
A1_D3_RUN7 = c(1110, 346, 919, 316, 652, 315)

# Manna-Whitney test
wilcox.test(A1_D3_old, A1_D3_RUN7)
# W = 13, p-value = 0.9143

A1_D6_old = c(271, 917)
A1_D6_RUN7 = c(1016, 219, 473, 426, 247)

# Manna-Whitney test
wilcox.test(A1_D6_old, A1_D6_RUN7)
# W = 6, p-value = 0.8571