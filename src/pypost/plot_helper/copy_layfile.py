import sys
from os import path, listdir
sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )


if __name__ == '__main__':
    ## Input:
    ## dir_in: directory for all the lay files
    ## filename_data: a dictionary with key as OLD data file, and value as NEW data file
    ## replaced:    if replaced orignial lay files
    dir_in = '/usr/local/home/yl8bc/Nic/PH/post/HLB_interp/'
    filename_data = {}
    filename_data['AveAcousticStat_timeave00137000-00217000_2D.plt'] = 'AveAcousticStat_timeave00137000-00337000_2D.plt'
    replaced = False
    ## Output:
    ## dir_out: output directory
    ## export_png: set True to export png to dir_out/
    dir_out = dir_in
    export_png = True
    ## NO MORE INPUT BELOW ##

    ## Get filenames under directory
    filenames_dir = listdir(dir_in)
    filenames_out_list = []
    ## Process lay files
    for filename_dir in filenames_dir:
        if '.lay' in filename_dir:
            f = open(dir_in+filename_dir, 'r')
            lay_old = f.read()
            f.close()
            lay_new = lay_old
            for data_old, data_new in filename_data.iteritems():
                lay_new = lay_new.replace(data_old, data_new)
            if replaced:
                filename_new = dir_out+filename_dir
            else:
                filename_new = dir_out+filename_dir[:-4]+'_new'+filename_dir[-4:]
            filenames_out_list.append(filename_new)
            f_new = open(filename_new, 'w')
            f_new.write(lay_new)
            f_new.close()

    ## generate .png files
    if export_png:
        import tecplot as tp
        import logging
        logging.basicConfig(level=logging.INFO)
        for filename in filenames_out_list:
            dset = tp.load_layout(filename)
            tp.export.save_png(filename.replace('.lay','.png'))

