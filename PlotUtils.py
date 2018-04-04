from __future__ import division, print_function
import os

try:
    from matplotlib.colors import hsv_to_rgb, LinearSegmentedColormap
    from matplotlib import cm
    from matplotlib import pyplot as mpl
    import pandas as pd
except:
    import tempfile
    import atexit
    import shutil

    mpldir = tempfile.mkdtemp()
    atexit.register(shutil.rmtree, mpldir)  # rm directory on succ exit

    os.environ['MPLCONFIGDIR'] = mpldir

    from matplotlib.colors import hsv_to_rgb, LinearSegmentedColormap
from scipy.stats import gaussian_kde
from numpy import array,  linspace, isfinite, median, exp, Inf
from itertools import repeat
import numpy as np
import subprocess
from collections import Counter

from os import path

ISH_ROT_4 = hsv_to_rgb(array(
    [[[(0.65+offset)%1, 0.00, 1.00],
      [(0.65+offset)%1, 0.53, 1.00],
      [(0.65+offset)%1, 0.53, 0.38],]
     for offset in linspace(0, 1, 4, endpoint=False)
    ]))
ISH_ROT_5 = hsv_to_rgb(array(
    [[[(0.65+offset)%1, 0.00, 1.00],
      [(0.65+offset)%1, 0.53, 1.00],
      [(0.65+offset)%1, 0.53, 0.38],]
     for offset in linspace(0, 1, 5, endpoint=False)
    ]))
ISH_ROT_6 = hsv_to_rgb(array(
    [[[(0.65+offset)%1, 0.00, 1.00],
      [(0.65+offset)%1, 0.53, 1.00],
      [(0.65+offset)%1, 0.53, 0.38],]
     for offset in linspace(0, 1, 6, endpoint=False)
    ]))

ISH_CMS_4 = []
ISH_CMS_5 = []
ISH_CMS_6 = []

for CMS, ROT in [(ISH_CMS_4, ISH_ROT_4),
                 (ISH_CMS_5, ISH_ROT_5),
                 (ISH_CMS_6, ISH_ROT_6)]:
    for I, ARR in enumerate(ROT):
        CMS.append(
            LinearSegmentedColormap('ish{}'.format(I),
                                    dict(red=((0.0, ARR[0, 0], ARR[0, 0]),
                                              (0.7, ARR[1, 0], ARR[1, 0]),
                                              (1.0, ARR[2, 0], ARR[2, 0])),
                                         green=((0.0, ARR[0, 1], ARR[0, 1]),
                                                (0.7, ARR[1, 1], ARR[1, 1]),
                                                (1.0, ARR[2, 1], ARR[2, 1])),
                                         blue=((0.0, ARR[0, 2], ARR[0, 2]),
                                               (0.7, ARR[1, 2], ARR[1, 2]),
                                               (1.0, ARR[2, 2], ARR[2, 2])),
                                        )))


ISH = LinearSegmentedColormap('ish',
                              dict(red=((0, 1, 1),
                                        (.7, 120/255, 120/255),
                                        (1, 46/255, 46/255)),
                                   green=((0, 1, 1),
                                          (.7, 129/255, 129/255),
                                          (1, 46/255, 46/255)),
                                   blue=((0, 244/255, 244/255),
                                         (.7, 1, 1),
                                         (1, 98/255, 98/255))))




def scatter_heat(x, y, **kwargs):
    kwargs['s'] = kwargs.get('s', 10)
    kwargs['edgecolors'] = kwargs.get('edgecolors', 'none')
    kwargs['cmap'] = kwargs.get('cmap', cm.jet)

    dropna = kwargs.pop('dropna', False)
    if dropna:
        ix = np.isfinite(x) & np.isfinite(y)
        x = x[ix]
        y = y[ix]

    if 'density' not in kwargs:
        estimator = gaussian_kde([x, y])
        density = estimator.evaluate([x, y])
    else:
        density = kwargs.pop('density')
    min_density = kwargs.pop('min_density', median(density))
    max_density = kwargs.pop('max_density', Inf)

    normdensity = exp(density.clip(min_density, max_density))
    xlim = kwargs.pop('xlim', (min(x), max(x)))
    ylim = kwargs.pop('ylim', (min(y), max(y)))
    retval = mpl.scatter(x, y, c=normdensity, **kwargs)
    ax = mpl.gca()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return retval, density

kwargs_single_column = dict(
    draw_box=True,
    draw_row_labels=False,
    max_width=221,
    total_width=200,
    vspacer=0,
    split_columns=True,
    col_sep='_sl',
    box_height=60,
    convert=True,
)

kwargs_heatmap = dict(
    draw_box=True,
    draw_row_labels=True,
    draw_name=True,
    total_width=90,
    vspacer=0,
    split_columns=True,
    col_sep='_sl',
    box_height=20,
    make_hyperlinks=True,
    convert=True,
)

kwargs_ase_heatmap = kwargs_heatmap.copy()
kwargs_ase_heatmap['cmap'] = cm.RdBu
kwargs_ase_heatmap['cmap'].set_bad((.5, .5, .5))
kwargs_ase_heatmap['norm_rows_by'] = 'center0pre'

kwargs_expr_heatmap = kwargs_heatmap.copy()
kwargs_expr_heatmap['norm_rows_by'] = 'maxall'

def svg_heatmap(data, filename, row_labels=None, box_size=4,
                index=None,
                all_indices=None, all_colnames=None, internal_datanames=None,
                cmap=ISH, norm_rows_by=None, draw_row_labels=False,
                color_row_labels=False,
                col_sep='', box_height=None, total_width=None,
                draw_box=False, draw_name=False, data_names=None,
                make_hyperlinks = False,
                progress_bar = False,
                max_width=np.inf,
                x_min=10, y_min=10,
                spacers=None,
                convert=False,
                squeeze_rows=None,
                cmap_by_prefix=None,
                draw_average=False,
                draw_average_only=False,
                average_scale=1,
                split_columns=False,
                vspacer=30,
                hatch_nan=True, hatch_size=20,
                figure_title=None,
                nan_replace=None,
                first_col='', last_col=''):
    """
    Draw heatmap as an SVG file stored in filename

    *data* can be either a 2D array-like type (list of lists, numpy array,
    pandas DataFrame, etc), or a tuple of 2D array-likes, in which case a
    separator will be added between each one in the output

    *cmap* is a matplotlib-like colormap (i.e. a callable that expects floats
    in the range 0.0-1.0.), or an iterable of the same length as the tuple
    *data* containing colormaps

    *row_labels* can be supplied, otherwise they will detected from the first
    item in *data*, if available, and if not they will be blank.

    If *total_width* is supplied, width of each dataset in *data* will be
    scaled to that constant. If *box_height* is supplied, the height of each
    row will be *box_height*, otherwise it will be equal to the width of each
    element. If neither are supplied, elements will be squares equal to
    *box_size*. IT IS STRONGLY RECOMMENDED that if if supplying *total_width*,
    *box_height* also be specified, but this is not enforced.

    *draw_row_labels*, if True, will label the rows on the right hand side. As
    of 2013/09/03, this won't scale the SVG properly, so including the
    resulting file in an html element won't display properly.

    *spacers* is the distance between adjacent datasets.  Can either be a
    number, in which case it will apply to all datasets, or an interable for
    different distances. If the iterable is shorter than the number of
    datasets, the last value will be repeated.

    """
    import svgwrite as svg
    try:
        import pandas as pd
        has_pandas = True
    except:
        has_pandas = False
        assert all_indices
        assert all_colnames

    if not isinstance(data, tuple):
        data = (data,)

    if not isinstance(norm_rows_by, tuple):
        norm_rows_by = repeat(norm_rows_by)


    old_data = data
    colname_tuple = repeat(None)
    if split_columns and has_pandas:
        from Utils import sel_startswith
        data = []
        new_normers = []
        new_cmaps = []
        if isinstance(cmap, tuple):
            cmaps = cmap
        else:
            cmaps = repeat(cmap)
        for dataset, normer, c_cmap in zip(old_data, norm_rows_by, cmaps):
            if dataset is None:
                data.append(dataset)
                new_normers.append(normer)
                new_cmaps.append(c_cmap)
                continue

            if not isinstance(dataset, pd.DataFrame):
                dataset = pd.DataFrame(dataset).T
            colnames = list(sorted(
                {col.split(col_sep)[0] for col in dataset.columns}))
            data.extend(
                dataset.select(**sel_startswith(colname)) for colname in colnames
            )
            new_normers.extend(normer for colname in colnames)
            new_cmaps.extend(c_cmap for colname in colnames)
        data = tuple(data)
        norm_rows_by = tuple(new_normers)
        cmap = tuple(new_cmaps)
    elif split_columns and all_colnames:
        colnames = list(sorted(
            {col.split(col_sep)[0] for col in all_colnames}))
        colname = colnames[0]
        data = tuple([
            data[:, array([c.startswith(colname) for c in internal_datanames])]
            for colname in colnames
        ])
        colname_tuple = tuple(
            [c for c in all_colnames if c.startswith(colname)]
            for colname in colnames
        )
    elif not split_columns and all_colnames:
        colname_tuple = tuple(
            [c for c in all_colnames if c.startswith(dataname)]
            for dataname in internal_datanames
        )

    rows, cols = np.shape([ds for ds in data if ds is not None][0])
    if index is not None:
        rows = len(index)
    if box_height is None:
        box_height = box_size

    if row_labels is None:
        if index is not None:
            row_labels = list(index)
        elif hasattr(data[0], 'index'):
            row_labels = list(data[0].index)
        else:
            row_labels = ['' for row in range(rows)]

    if total_width is not None and max_width is not np.inf:
        boxes_per_row = max_width // (1.1 * total_width)
        if ((boxes_per_row + 1) * 1.1 * total_width - .1 * total_width
            < max_width):
            boxes_per_row += 1

        num_plotted_rows = np.ceil(len(data) / boxes_per_row
                                   + (draw_average or draw_average_only))
        if figure_title is None:
            fig_title_height = 0
        elif isinstance(figure_title, tuple):
            fig_title_height = len(figure_title)
        else:
            fig_title_height = 1
        dwg = svg.Drawing(filename,
                          size=(max_width + 2 * x_min + 200 * draw_row_labels,
                                2 * y_min
                                + (num_plotted_rows
                                   * (rows)
                                   * box_height)
                                + 80 * (fig_title_height)
                                + 80 * draw_name
                                + (num_plotted_rows - 1) * vspacer),
                         )
    elif total_width is not None:
        width = len(data) * total_width * 1.1 - .1 * total_width
        height = rows * box_height
        max_row_label_len = max(len(str(i)) for i in row_labels)
        dwg = svg.Drawing(filename,
                          size=(width + 2 * x_min + 20 * draw_row_labels *
                                max_row_label_len,
                                height + 2 * y_min + 80 * draw_name
                                + (80 * (figure_title is not None)))
                         )
    else:
        dwg = svg.Drawing(filename)
    dwg.add(svg.base.Title(path.basename(filename)))

    pat = dwg.pattern(id='hatch', insert=(0, 0), size=(hatch_size, hatch_size),
                      patternUnits='userSpaceOnUse')
    g = pat.add(dwg.g(style="fill:none; stroke:#B0B0B0; stroke-width:1"))
    g.add(dwg.path(('M0,0', 'l{hatch},{hatch}'.format(hatch=hatch_size))))
    g.add(dwg.path(('M{hatch2},0 l{hatch2},{hatch2}'.format(hatch2=hatch_size/2).split())))
    g.add(dwg.path(('M0,{hatch2} l{hatch2},{hatch2}'.format(hatch2=hatch_size/2).split())))

    dwg.add(pat)


    if box_height is None:
        box_height = box_size

    if not hasattr(cmap, "__len__"):
        cmap = [cmap for frame in data]

    if data_names is None:
        data_names = ["" for frame in data]

    if len(cmap) != len(data):
        raise ValueError("cmap and data should be the same length ({} vs {})"
                        .format(len(cmap), len(data)))

    if not hasattr(spacers, "__len__"):
        spacers = [spacers]
    else:
        spacers = list(spacers)
    while len(spacers) < len(data):
        spacers.append(spacers[-1])

    if ((isinstance(norm_rows_by, repeat)
         and isinstance(next(norm_rows_by), str)
         and next(norm_rows_by).startswith('center0all'))
        or (not isinstance(norm_rows_by, repeat)
            and isinstance(norm_rows_by[0], str)
            and np.any([i.startswith('center0all') for i in norm_rows_by]))):
        all_data = pd.concat(data, axis=1)

    if squeeze_rows is not None:
        data = [
            pd.DataFrame(d.apply(squeeze_rows, axis=1),
                         columns=[path.commonprefix(list(d.columns))])
            for d in data
        ]

    x_start = x_min
    y_start = y_min
    y_diff = 0
    iterator = zip(data, cmap, data_names, norm_rows_by, spacers,
                   colname_tuple)
    if figure_title:
        if isinstance(figure_title, tuple):
            font_size = '3em'
            for title_line in figure_title:
                dwg.add(dwg.text(title_line, (x_start, y_start+75,),
                                 style="font-size:{};font-family:sans-serif".format(font_size)))
                y_start += 80
                font_size = '1.5em'

        else:
            dwg.add(dwg.text(figure_title, (x_start, y_start+75,),
                             style="font-size:3em;font-family:sans-serif"))
            y_start += 80
    if progress_bar:
        from progressbar import ProgressBar
        pbar = ProgressBar(maxval=len(data)*rows).start()
        pbar_val = 0

    for frame, c_cmap, name, normer, spacer, colnames in iterator:
        if frame is None:
            dwg.add(dwg.text(normer, (x_start, y_start + box_height/2)))
            if total_width is not None:
                if spacer is None:
                    x_start += total_width * 1.1
                else:
                    x_start += total_width + spacer
            else:
                if spacer is None:
                    x_start += box_size
                else:
                    x_start += spacer
            if x_start > max_width:
                x_start = x_min
                y_start += box_height + vspacer
            continue
        if has_pandas:
            frame = pd.DataFrame(frame)
        if index is not None:
            if has_pandas:
                frame = frame.ix[index]
            else:
                setix = set(index)
                #selector = [i for i, name in enumerate(all_indices) if name in setix]
                #frame = frame[selector, :]
        if normer is None:
            norm_data = array(frame.copy())
        elif normer is 'mean':
            if has_pandas:
                norm_data = array(frame.divide(frame.dropna(axis=1, how='all').mean(axis=1)+10, axis=0))
            else:
                norm_data = frame / (frame[:,isfinite(frame[0,:])].mean(axis=1) + 10).reshape((rows, 1))
        elif normer == 'max':
            if has_pandas:
                norm_data = array(frame.divide(frame.dropna(axis=1, how='all').max(axis=1)+10, axis=0))
            else:
                norm_data = frame / (frame[:,isfinite(frame[0,:])].max(axis=1) + 10).reshape((rows, 1))
        elif normer == 'maxall':
            if has_pandas:
                maxall = frame.max(axis=1)
                assert len(data) == len(new_normers)
                for old_frame, norm_type in zip(data, new_normers):
                    if norm_type != 'maxall': continue
                    if old_frame is not None:
                        old_frame = old_frame.max(axis=1).ix[index
                                                             if index is not None
                                                             else old_frame.index]
                        maxall = maxall.where(maxall > old_frame, old_frame)
                norm_data = array(frame.divide(maxall + 10, axis=0))
            else:
                norm_data = frame / (old_data[:, isfinite(old_data[0, :])]
                                     .max(axis=1) + 10).reshape((rows, 1))
        elif normer == 'fullvar':
            norm_data = frame.subtract(frame
                                       .dropna(axis=1, how='all')
                                       .min(axis=1)-1e-6,
                                       axis=0)
            norm_data = array(norm_data.divide(norm_data
                                               .dropna(axis=1, how='all')
                                               .max(axis=1),
                                               axis=0))
        elif normer == 'center0':
            norm_data = array(0.5 +
                         0.5 * frame.divide(frame.dropna(axis=1).abs().max(axis=1),
                                      axis=0)
                        )
        elif isinstance(normer, str) and normer.startswith('center0min'):
            min_norm = (
                frame.dropna(axis=1).abs() .max(axis=1).clip(float(normer[10:]), 1e6)
            )
            norm_data = array(0.5+
                              0.5 * frame.divide(min_norm, axis=0))

        elif isinstance(normer, str) and normer.startswith('center0allmin'):
            min_norm = (
                all_data.dropna(axis=1).abs() .max(axis=1).clip(float(normer[13:]), 1e6)
            )
            norm_data = array(0.5+
                              0.5 * frame.divide(min_norm, axis=0))

        elif normer == 'center0all':
            norm_data = array(0.5 +
                         0.5 *
                         frame.divide(all_data.dropna(how='all', axis=1).abs().max(axis=1),
                                      axis=0)
                        )
        elif normer == 'center0pre':
            norm_data = array(0.5 + 0.5 * frame)
        elif isinstance(normer, (int, float)):
            norm_data = array(frame / normer)
            normer = 'constant'
        elif index is not None and hasattr(normer, "ix"):
            norm_data = array(frame.divide(normer.ix[index], axis=0))
        elif hasattr(normer, "__len__") and len(normer) == rows:
            if has_pandas:
                norm_data = array(frame.divide(normer, axis=0))
            else:
                norm_data = array(frame / np.reshape(normer, (rows, -1)))


        elif hasattr(normer, "__len__"):
            print('\n'*5)
            print(len(normer), normer, normer=='max')
            print(frame.shape)
            raise TypeError("norm_rows_by should be the same shape "
                            "as the number of rows")
        else:
            norm_data = array(frame / normer)

        if not c_cmap or str(c_cmap).lower() == 'default':
            c_cmap = ISH

        new_rows, new_cols = np.shape(frame)
        if hasattr(frame, 'index'):
            col_labels = frame.columns
        elif colnames:
            col_labels = colnames
        else:
            col_labels = ['' for col in range(new_cols)]
        if new_rows != rows:
            raise ValueError("All input elements must have the same number of"
                             " rows (and same row meanings --unchecked)")

        if total_width is not None:
            box_size = total_width / float(new_cols)

        i = 0
        if not draw_average_only:
            for i in range(rows):
                if progress_bar:
                    pbar.update(pbar_val)
                    pbar_val += 1
                prefix = col_labels[0][:col_labels[0].find(col_sep)]
                if cmap_by_prefix:
                    c_cmap = cmap_by_prefix(prefix)
                for j in range(new_cols):
                    g = dwg.g()
                    val = frame.ix[i,j] if has_pandas else frame[i,j]
                    g.add(svg.base.Title("{}, {}: {:.2f}".format(row_labels[i],
                                                                 col_labels[j],
                                                                 val)))
                    hatch = not isfinite(norm_data[i,j])
                    if hatch and nan_replace is not None:
                        if isinstance(nan_replace, float):
                            norm_data[i, j] = nan_replace
                        else:
                            if normer.startswith('center0'):
                                norm_data[i, j] = 0.5
                            else:
                                norm_data[i, j] = 0.0
                    elif hatch:
                        n = 0
                        norm_data[i, j] = 0
                        left = j - 1
                        while left >= 0:
                            if isfinite(norm_data[i, left]):
                                norm_data[i, j] += norm_data[i, left]
                                n += 1
                                break
                            left -= 1
                        right = j + 1
                        while right  < norm_data.shape[1]:
                            if isfinite(norm_data[i, right]):
                                norm_data[i, j] += norm_data[i, right]
                                n+= 1
                                break
                            right += 1
                        if n == 0:
                            norm_data[i, j] = .5 if 'center' in normer else 0
                        else:
                            norm_data[i, j] /= n
                    g.add(dwg.rect((x_start + box_size*j, y_start + i*box_height),
                                   (box_size, box_height),
                                   style="fill:#{:02x}{:02x}{:02x}"
                                   .format(*[int(255*x) for x in
                                             c_cmap(norm_data[i, j])])))
                    dwg.add(g)
                    if hatch_nan and hatch:
                        g.add(dwg.rect((x_start + box_size*j,
                                        y_start + i*box_height),
                                       (box_size, box_height),
                                       style="fill:url(#hatch)"
                                      )
                             )
                    col_base = col_labels[j][:col_labels[j].find(col_sep)]
                    if col_base != prefix:
                        prefix = col_base
                        if cmap_by_prefix:
                            c_cmap = cmap_by_prefix(prefix)
                        g.add(dwg.line((x_start + box_size * j,
                                        y_start + i * box_height),
                                       (x_start + box_size * j,
                                        y_start + (i + 1) * box_height),
                                       style="stroke-width:{}; stroke:#000000"
                                       .format(.1 * box_size)))
        else:
            for j in range(new_cols):
                hatch = not isfinite(norm_data[0, j])
                if hatch:
                    n = 0
                    norm_data[:, j] = 0
                    if j > 0 and isfinite(norm_data[0,j-1]):
                        norm_data[:, j] += norm_data[:, j-1]
                        n += 1
                    if (j + 1 < norm_data.shape[1]
                        and isfinite(norm_data[0, j+1])):
                        norm_data[:, j] += norm_data[:, j+1]
                        n += 1
                    norm_data[:, j] /= n
        dwg.add(dwg.text(first_col, (x_start,
                                     y_start + (i + 1) * box_height)))
        dwg.add(dwg.text(last_col, (x_start + (new_cols - 1) * box_size,
                                    y_start + (i + 1) * box_height)))
        if draw_box and not draw_average_only:
            dwg.add(dwg.rect((x_start, y_start + 0),
                             (new_cols*box_size, rows*box_height),
                             style="stroke-width:1; "
                             "stroke:#000000; fill:none"))
        if draw_average or draw_average_only:
            avg_frame = norm_data.mean(axis=0)
            for j in range(new_cols):
                col_base = col_labels[j][:col_labels[j].find(col_sep)]
                prefix = col_base
                if cmap_by_prefix:
                    c_cmap = cmap_by_prefix(prefix)
                g = dwg.g()
                g.add(svg.base.Title("Average, {}: {:.2f}".format(col_labels[j],
                                                                  avg_frame[j])))
                g.add(dwg.rect((x_start + box_size*j,
                                y_start + (i+(not draw_average_only))*box_height),
                               (box_size, box_height),
                               style="fill:#{:02x}{:02x}{:02x}"
                               .format(*[int(255*x) for x in
                                         c_cmap(average_scale*avg_frame[j])])))
                if not isfinite(norm_data[0, j]) and hatch_nan:
                    g.add(dwg.rect((x_start + box_size*j,
                                    y_start + (i+(not draw_average_only))*box_height),
                                   (box_size, box_height),
                                   style="fill:url(#hatch)"
                                  )
                         )

                dwg.add(g)
            dwg.add(dwg.rect((x_start,
                              y_start + (i+(not draw_average_only))*box_height),
                             (new_cols*box_size, 1*box_height),
                             style="stroke-width:1; stroke:#000000; fill:none"
                            ))


        if draw_name:
            if name == "" and split_columns:
                name = col_base
            xpos = x_start + box_size * new_cols / 2.0
            text = dwg.text('',
                             (xpos,
                              y_start
                              + box_height * (rows) * (1-draw_average_only)
                              + box_height * (draw_average or draw_average_only)
                              + 13),
                            style="text-anchor: middle;font-family:sans-serif;")
            text.add(dwg.tspan("", dy=["-1.5em"]))
            for line in name.split('_'):
                text.add(dwg.tspan(line,
                                   dy=["1.5em"],
                                   x=[xpos],
                                   style="text-anchor: middle;",
                                   ))
            dwg.add(text)

        if total_width is not None:
            if spacer is None:
                x_start += total_width * 1.1
            else:
                x_start += total_width + spacer
        else:
            if spacer is None:
                x_start += new_cols * box_size + box_size
            else:
                x_start += new_cols * box_size + spacer

        #y_diff = new_rows * box_height + vspacer
        if x_start + total_width >= max_width:
            x_start = x_min
            y_start += new_rows*box_height*(not draw_average_only) + vspacer
            y_start += box_height*(draw_average_only or draw_average)

    if draw_row_labels and isinstance(row_labels[0], tuple):
        lwidths = Counter()
        for r in row_labels:
            for i, l in enumerate(r):
                lwidths[i] = max(lwidths[i], len(str(l)))
        cum_len = 0
        for i in range(len(lwidths)):
            old_width = lwidths[i]
            lwidths[i] += cum_len
            cum_len += old_width

    if draw_row_labels and not draw_average_only:
        for i in range(rows):
            if color_row_labels:
                style = "font-family:sans-serif; font-size: {size}; fill: {color};".format(
                    size=box_height,
                    color='red' if row_labels[i] in color_row_labels else 'black',
                )
            else:
                style = "font-family:sans-serif; font-size: {}".format(box_height)
            if isinstance(row_labels[i], tuple):
                labeltext = dwg.g()
                for lnum, ltext in enumerate(row_labels[i]):
                    labeltext.add(dwg.text(ltext,
                                           (x_start + lwidths[lnum-1] * 10 + lnum * 50,
                                            y_start + i * box_height + box_height),
                                           style=style,
                                          ))
            else:
                labeltext = (dwg.text(row_labels[i],
                                      (x_start, y_start + i*box_height+box_height),
                                      style=style,
                                     ))
            if make_hyperlinks:
                if make_hyperlinks is True:
                    link = dwg.a('http://insitu.fruitfly.org/cgi-bin/ex/report.pl?ftype={}&ftext={}'
                                 .format(2 if (isinstance(row_labels[i], str)
                                               and
                                               (row_labels[i].startswith('FBgn'))
                                              )
                                         else 1,
                                         row_labels[i]),
                                 target='_replace',
                                )
                else:
                    link = dwg.a(make_hyperlinks.format(frame.index[i]))
                link.add(labeltext)
                dwg.add(link)
            else:
                dwg.add(labeltext)
    if progress_bar:
        pbar.finish()
    dwg.saveas(filename)
    if convert:
        cmd = [
            'convert',
            filename,
            '-units', 'PixelsPerInch',
            '+antialias',
            '-density', '600',
            '-background', 'none',
            '-transparent', 'white',
            filename.replace('svg', 'png'),
        ]
        subprocess.Popen(cmd)



def cmap_by_prefix(prefix):
    cms = dict(
        WT = ISH_CMS_5[0],
        bcd = ISH_CMS_5[1],
        zld = ISH_CMS_5[2],
        G20 = ISH_CMS_5[3],
        hb = ISH_CMS_5[4],
    )
    for p in cms:
        if prefix.startswith(p):
            return cms[p]
    return ISH

def minimize_ink(axes):
    axes.spines['right'].set_color('none')
    axes.spines['top'].set_color('none')
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')


def plot_coordinates(chroms, coords, **kwargs):
    offsets = pd.Series({
        '2L': 23011544,
        '2R': 0,
        '3L': 24543557,
        '3R': 0,
        '4': 0,
        'X': 0
    })
    chrom_lens = {
        "YHet"                       :347038,
        "2RHet"                      :3288761,
        "2LHet"                     :368872,
        "3LHet"                      :2555491,
        "3RHet"                      :2517507,
        "U"                          :10049037,
        "XHet"                       :204112,
        "dmel_mitochondrion_genome"  :19517,
        "2L"                         :23011544,
        "X"                          :22422827,
        "3L"                         :24543557,
        "4"                          :1351857,
        "2R"                         :21146708,
        "3R"                         :27905053,
        "Uextra"                     :29004656,
    }

    chrom_y_coords = pd.Series({
        '2L': 1,
        '2R': 1,
        '3L': 2,
        '3R': 2,
        '4': 3,
        'X': 4,
    })

    for chrom in chrom_y_coords.index:
        mpl.hlines(chrom_y_coords[chrom], chrom_lens[chrom] - offsets[chrom],
                   -offsets[chrom])
    return mpl.scatter(
        coords - pd.Series(index=chroms.index, data=offsets[chroms].as_matrix()),
        chrom_y_coords[chroms],
        **kwargs
    )
