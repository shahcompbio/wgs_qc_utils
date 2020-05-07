breakpoints = pd.read_csv(breakpoints)[["chromosome_1", "chromosome_2", "position_1", "position_2"]]
breakpoints = breakpoints.astype({"chromosome_1": str, "chromosome_2": str})

breakpoints = pd.DataFrame({"chr": breakpoints["chromosome_1"].append(breakpoints["chromosome_2"]),
                            "pos": breakpoints["position_1"].append(breakpoints["position_2"])})
