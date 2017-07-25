#!/usr/bin/env python3
import argparse

from root_pandas import read_root


tree_key = 'genResults'
columns = ['m12sq', 'm13sq', 'm24sq', 'm34sq', 'm134sq']


def prepare_data(input_fn, output_fn_d0, output_fn_d0bar):
    df = read_root(input_fn, key=tree_key)

    df_d0 = df.query('isAntiD0 == 0')[columns]
    df_d0_bar = df.query('isAntiD0 == 1')[columns]

    df_d0.to_csv(output_fn_d0, header=False, index=False, sep=' ')
    df_d0_bar.to_csv(output_fn_d0bar, header=False, index=False, sep=' ')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-fn', required=True)
    parser.add_argument('--output-fn-D0', required=True)
    parser.add_argument('--output-fn-D0bar', required=True)
    args = parser.parse_args()
    prepare_data(args.input_fn, args.output_fn_D0, args.output_fn_D0bar)
