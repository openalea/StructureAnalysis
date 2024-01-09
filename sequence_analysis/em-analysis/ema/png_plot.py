"""
.. module:: png_plot
   :platform: Unix, Windows
   :synopsis: Visualisation et sauvegarde de scanpaths

"""

from config import COLOR_PALETTE
from config import PNG_TEXTS_PATH
import numpy as np
import os
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont


def seaborn_color_palette_to_pil():
    """ Converts a list of sns rgb colors to pil rgb colors

    :return: list of str ['rgb(10%, 10%, 10%)', ..., ]
    """
    pil_color_palette = []

    for col in xrange(0, np.shape(COLOR_PALETTE)[0]):
        red = int(round(COLOR_PALETTE[col][0] * 100))
        green = int(round(COLOR_PALETTE[col][1] * 100))
        blue = int(round(COLOR_PALETTE[col][2] * 100))
        pil_color_palette.append('rgb(' + str(red) + '%, ' + str(green) + '%, ' + str(blue) + '%)')

    return pil_color_palette


def seq_id_dic(df):
    """
    Return a dictionary {seq_num: (subj_name, text)}
    """
    d = {}
    s = 0

    for last_fixation_index in df[df.ISLAST == 1].index:  # browse all seq
        d[s] = (df.at[last_fixation_index, 'SUBJ_NAME'], df.at[last_fixation_index, 'TEXT'])

    return d


def plot_scanpath(df, subj, text, output_file="", display=True):
    """Display scanpath without state restoration.

    Args:
        df: pandas dataframe
        subj: subject name (e.g. s01)
        text: text name (e.g. aide_refugies-f1)
        output_file: if given then img is saved
        display: display img or not. Do it by default
    """
    assert df[(df.SUBJ_NAME == subj) & (df.TEXT == text)].ISLAST.sum() == 1
    assert df[(df.SUBJ_NAME == subj) & (df.TEXT == text)].ISFIRST.sum() == 1

    img = Image.open(os.path.join(PNG_TEXTS_PATH, text) + '.png')
    draw = ImageDraw.Draw(img)

    radius_scale = 0.1  # radius scaling factor
    c = 'black'

    for fixation_index in df[(df.SUBJ_NAME == subj) & (df.TEXT == text)].index:  # browse all seq
        if df.at[fixation_index, 'ISFIRST'] != 1:
            previous_x = x
            previous_y = y
        x = int(df.at[fixation_index, 'X'])
        y = int(df.at[fixation_index, 'Y'])
        fdur = df.at[fixation_index, 'FDUR']
        r = int(radius_scale * fdur)
        draw.ellipse([x - r, y - r, x + r, y + r], outline='black')
        if df.at[fixation_index, 'ISFIRST'] != 1:
            draw.line([previous_x, previous_y, x, y], fill='black')

    if output_file != "":
        img.save(output_file)
    if display:
        img.show()

    return img


def plot_state_scanpath(df, subj, text, output_file="", display=True):
    """Display scanpath with state restoration.

    Args:
        df: pandas dataframe
        subj: subject name (e.g. s01)
        text: text name (e.g. aide_refugies-f1)
        output_file: if given then img is saved
        display: display img or not. Do it by default
    """
    assert df[(df.SUBJ_NAME == subj) & (df.TEXT == text)].ISLAST.sum() == 1
    assert df[(df.SUBJ_NAME == subj) & (df.TEXT == text)].ISFIRST.sum() == 1
    assert 'STATES' in df.columns

    img = Image.open(os.path.join(PNG_TEXTS_PATH, text) + '.png')
    draw = ImageDraw.Draw(img)
    color_list = seaborn_color_palette_to_pil()
    radius_scale = 0.1  # radius scaling factor

    for fixation_index in df[(df.SUBJ_NAME == subj) & (df.TEXT == text)].index:  # browse all seq
        if df.at[fixation_index, 'ISFIRST'] != 1:
            previous_x = x
            previous_y = y
            previous_state = state
        x = int(df.at[fixation_index, 'X'])
        y = int(df.at[fixation_index, 'Y'])
        state = int(df.at[fixation_index, 'STATES'])
        fdur = df.at[fixation_index, 'FDUR']
        r = int(radius_scale * fdur)
        draw.ellipse([x - r, y - r, x + r, y + r], outline=color_list[state])
        if df.at[fixation_index, 'ISFIRST'] != 1:
            draw.line([previous_x, previous_y, x, y], fill=color_list[previous_state])

    if output_file != "":
        img.save(output_file)
    if display:
        img.show()
    return img


def plot_col_scanpath(df, col, subj, text, output_file="", display=True):
    """Display scanpath with col for each fixation (e.g. 'READMODE' or ['CINC', 'WINC']).

    Args:
        df: pandas dataframe
        col: a column or a set of columns
        subj: subject name (e.g. s01)
        text: text name (e.g. aide_refugies-f1)
        output_file: if given then img is saved, otherwise it is displayed
        display: display img or not. Do it by default
    """
    assert df[(df.SUBJ_NAME == subj) & (df.TEXT == text)].ISLAST.sum() == 1
    assert df[(df.SUBJ_NAME == subj) & (df.TEXT == text)].ISFIRST.sum() == 1
    assert type(col) == str or type(col) == list
    if type(col) == str:
        assert col in df.columns
    elif type(col) == list:
        for e in col:
            assert type(e) == str
            assert e in df.columns

    img = Image.open(os.path.join(PNG_TEXTS_PATH, text) + '.png')
    draw = ImageDraw.Draw(img)

    radius_scale = 0.1  # radius scaling factor
    c = 'black'
    font = ImageFont.load_default().font

    for fixation_index in df[(df.SUBJ_NAME == subj) & (df.TEXT == text)].index:  # browse all seq
        if df.at[fixation_index, 'ISFIRST'] != 1:
            previous_x = x
            previous_y = y
        x = int(df.at[fixation_index, 'X'])
        y = int(df.at[fixation_index, 'Y'])
        fdur = df.at[fixation_index, 'FDUR']
        r = int(radius_scale * fdur)
        draw.ellipse([x - r, y - r, x + r, y + r], outline='black')
        if type(col) == str:
            draw.text((x, y), str(df.at[fixation_index, col]), 'black', font=font)
        elif type(col) == list:
            shift = 0
            for e in col:
                draw.text((x, y + shift), str(df.at[fixation_index, e]), 'black', font=font)
                shift += 10
        if df.at[fixation_index, 'ISFIRST'] != 1:
            draw.line([previous_x, previous_y, x, y], fill='black')

    if output_file != "":
        img.save(output_file)
    if display:
        img.show()
    return img


def plot_state_and_col_scanpath(df, col, subj, text, output_file="", display=True):
    """Display scanpath with state restoration and col for each fixation (e.g. 'READMODE' or ['CINC', 'WINC']).

        Args:
        df: pandas dataframe
        col: a column or a set of columns
        subj: subject name (e.g. s01)
        text: text name (e.g. aide_refugies-f1)
        output_file: if given then img is saved, otherwise it is displayed
        display: display img or not. Do it by default
    """

    assert df[(df.SUBJ_NAME == subj) & (df.TEXT == text)].ISLAST.sum() == 1
    assert df[(df.SUBJ_NAME == subj) & (df.TEXT == text)].ISFIRST.sum() == 1
    assert 'STATES' in df.columns
    assert type(col) == str or type(col) == list
    if type(col) == str:
        assert col in df.columns
    elif type(col) == list:
        for e in col:
            assert type(e) == str
            assert e in df.columns

    img = Image.open(os.path.join(PNG_TEXTS_PATH, text) + '.png')
    draw = ImageDraw.Draw(img)
    color_list = seaborn_color_palette_to_pil()
    radius_scale = 0.1  # radius scaling factor
    font = ImageFont.load_default().font

    for fixation_index in df[(df.SUBJ_NAME == subj) & (df.TEXT == text)].index:  # browse all seq
        if df.at[fixation_index, 'ISFIRST'] != 1:
            previous_x = x
            previous_y = y
            previous_state = state
        x = int(df.at[fixation_index, 'X'])
        y = int(df.at[fixation_index, 'Y'])
        state = int(df.at[fixation_index, 'STATES'])
        fdur = df.at[fixation_index, 'FDUR']
        r = int(radius_scale * fdur)
        draw.ellipse([x - r, y - r, x + r, y + r], outline=color_list[state])
        if type(col) == str:
            draw.text((x, y), str(df.at[fixation_index, col]), 'black', font=font)
        elif type(col) == list:
            shift = 0
            for e in col:
                draw.text((x, y + shift), str(df.at[fixation_index, e]), 'black', font=font)
                shift += 10
        if df.at[fixation_index, 'ISFIRST'] != 1:
            draw.line([previous_x, previous_y, x, y], fill=color_list[previous_state])

    if output_file != "":
        img.save(output_file)
    if display:
        img.show()
    return img
