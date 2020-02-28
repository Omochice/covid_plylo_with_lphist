import numpy as np
from bresenham import bresenham
import itertools


def culc_local_binary_pattern(filled_boolean_list):
    # [T, F, T, T, F, T, F, T, F]みたいなリストを受け取る
    pattern_number = 0
    for i, part in enumerate(filled_boolean_list):
        if part:
            pattern_number += 2**(8 - i)
    return pattern_number


def generate_mini_windows(filled_pixels):
    previous_top_left_points = set([])
    for filled_pixel in filled_pixels:
        for index in range(9):  # 左上to右下
            top_left = (filled_pixel[0] - (index % 3),
                        filled_pixel[1] + (index // 3))
            if top_left not in previous_top_left_points:  #左上が決まれば小窓は一意に決まる
                previous_top_left_points |= set([top_left])  # 使った左上を保存
                yield generate_mini_window(top_left)


def generate_mini_window(top_left: tuple):
    window = [(top_left[0] + dx, top_left[1] - dy)
              for dy, dx in itertools.product(range(3), repeat=2)]
    # 以下と等価
    # row1 = [(top_left[0] + dx, top_left[1]) for dx in range(3)]
    # row2 = [(top_left[0] + dx, top_left[1] - 1) for dx in range(3)]
    # row3 = [(top_left[0] + dx, top_left[1] - 2) for dx in range(3)]
    # window = row1 + row2 + row3
    return window


def make_lphist(x_coodinates: list, y_coodinates: list,
                normalize: bool = True) -> list:
    coodinates = [(x, y) for x, y in zip(x_coodinates, y_coodinates)]

    filled_pixels = set([])

    results = np.zeros(512, dtype=np.uint8)

    for first, last in zip(coodinates[:-1], coodinates[1:]):
        first_x = int(first[0])
        first_y = int(first[1])
        last_x = int(last[0])
        last_y = int(last[1])

        # setで和をとる
        filled_pixels |= set(bresenham(first_x, first_y, last_x, last_y))

    for mini_window in generate_mini_windows(filled_pixels):
        boolean_is_filled_list = [(pixel in filled_pixels)
                                  for pixel in mini_window]
        results[culc_local_binary_pattern(boolean_is_filled_list)] += 1

    if normalize:
        results = results / np.sum(results)
    return results


if __name__ == "__main__":
    x_coo = [1.0, 2.0, 3.547698135744552, 1.5826141681866683]
    y_coo = [-1.0, 0.0, 1.547698135744552, 3.5127821033024356]
    for index, value in enumerate(make_lphist(x_coo, y_coo)):
        if value > 0:
            print(index, value)

    # print(make_local_pattern_histgram(x_coo, y_coo))
