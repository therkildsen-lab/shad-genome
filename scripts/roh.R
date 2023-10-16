bcftools_roh_042123 %>%
    filter(`# ST` == "RG") %>%
    mutate(length = `[4]Position` - `[5]State (0:HW, 1:AZ)`)
