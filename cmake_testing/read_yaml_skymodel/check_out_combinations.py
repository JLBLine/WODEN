

comp_types = ['POINT', 'GAUSSIAN', 'SHAPELET']
source_configs = ['SOURCE', 'MULTI_SOURCE']
comp_configs = ['COMP', 'MULTI_COMP']
flux_types = ['POWER', 'CURVE', 'LIST']


# for comp_type in comp_types:
#     for source_config in source_configs:
#         for comp_config in comp_configs:
#             for flux_type in flux_types:
#
#                 print(f"{comp_type} {source_config} {comp_config} {flux_type}")


for comp_type in comp_types:
    for flux_type in flux_types:
        print(f"{comp_type} {flux_type}")
