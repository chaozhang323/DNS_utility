import numpy as np
import os

class wave:
    """A class of wave recognition """
    dirs_dict = {'pics':'recognized_pics/', 'save':'save/'}
    flnm_mark_dict = {'.npy':'bin_', '.png':'pic_'}

    def __init__(self, name, path_pics, path_result):
        self.name = name
        self.path_pics = path_pics
        self.path_result = path_result

    def update_parameters(self, **kwargs):
        self.C_threshold = kwargs['C_T']
        self.morph_ele = kwargs['morph_ele']
        self.C_cut = kwargs['C_cut']
        self.C_complet = kwargs['C_complete']
        self.C_border = kwargs['C_border']
        self.num_frame_skipped = kwargs['num_frame_skipped']
        self.C_corr = kwargs['C_corr']

    def gather_pics_filenames(self, format):
        flnms = os.listdir(self.path_pics)
        flnms_pics = []
        for flnm in flnms:
            if format in flnm:
                flnms_pics.append(flnm)
        return sorted(flnms_pics)

    def make_dir(self):
        for dir in self.dirs_dict.values():
            if not os.path.exists(self.path_result+dir):
                os.mkdir(self.path_result+dir)




w = wave('M8', '/usr/local/home/yl8bc/PH/post/acoustic_radiation/M8/pics/',  './acoustic_radiation/M8/')
w.gather_pics_filenames('.png')