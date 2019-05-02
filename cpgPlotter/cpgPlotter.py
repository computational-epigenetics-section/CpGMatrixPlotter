import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class CpGMatrixPlotter:
    """
    Main class used to plot CpG matrix data
    Example:
    >>>from cpgPlotter import CpgMatrixPlotter
    >>>import numpy as np
    >>>data = np.array([[1,1,1,0],
    >>>                 [1,1,0,0]])
    >>>locations = np.array([100, 110, 155, 190])
    >>>plotter = CpgMatrixPlotter()
    >>>plotter.plotCpgMatrix(data, locations)
    """

    def __init__(self, highlight_color="limegreen"):
        self.highlight_color = highlight_color

    @staticmethod
    def _buffer_spacings(spacings, radius):
        for i in range(len(spacings)):
            if i == 0:
                continue
            distance = spacings[i] - spacings[i - 1]
            if distance < radius:
                spacings[i] = spacings[i] + (radius - distance)

        return np.array(spacings)

    @staticmethod
    def _get_range_from_bin(bin_str: str, bin_size=100):
        chrom, loc = bin_str.split("_")
        loc = int(loc)
        start = loc - bin_size
        return start, loc

    @staticmethod
    def prep_clustered_data_frame(clustered_data_frame):
        cluster_labels = clustered_data_frame['class']
        input_labels = clustered_data_frame['input']

        clustered_data_frame['meth_mean'] = clustered_data_frame.drop(['input', 'class'], axis=1).mean(axis=1)
        clustered_data_frame = clustered_data_frame.sort_values(['meth_mean', 'class'], ascending=False)

        working_df = clustered_data_frame.drop(['class', 'input', 'meth_mean'], axis=1)

        cpgMatrix = np.array(working_df)
        cpgPositions = np.array([int(x) for x in working_df.columns])

        return cpgMatrix, cpgPositions, cluster_labels, input_labels

    @staticmethod
    def _get_color(cpg_value: int):
        if cpg_value == 1:
            return "black"
        elif cpg_value == 0:
            return "white"

        else:
            NotImplementedError("I cannot yet accept unknown values. But I will soon.")

    @staticmethod
    def _sort_matrix_by_methylation(cpgMatrix):
        df = pd.DataFrame(cpgMatrix)
        df['mean'] = df.apply(np.mean, axis=1)
        df = df.sort_values(['mean'], ascending=True)

        return np.array(df.drop(['mean'], axis=1))
    
    def _get_highlight_color(self, cpg_highlight: int):
        if cpg_highlight == 1:
            return self.highlight_color
        else:
            return "black"

    def plotCpGMatrix(self, cpgMatrix, cpgPositions, title=None, figsize=(6, 6), sort=False, highlights=None):
        """
        Main plotting function. Plot tanghulu plot of cpg matrix data
        :param cpgMatrix: np.array of CpG values. 1=methylated; 0=unmethylated
        :param cpgPositions: genomic positions of the CpGs. Used to calculate spacing
        :param title: Title to add to the plot, optional
        :param figsize: size of the plot to generate, default (8,8)
        :param sort: sort the reads my methylation value
        :return: matplotlib axes
        """
        fig, ax = plt.subplots(figsize=figsize)
        ax.axis("equal")
        v_steps = 1 / cpgMatrix.shape[0]
        v_spacings = np.arange(0, 1, v_steps)
        h_spacings = (cpgPositions - min(cpgPositions)) * 0.01
        ax.set_ylim(-.05, 1.05)
        ax.set_xlim(-.02, max(h_spacings) + 0.02)
        ax.set_xticks(h_spacings)
        ax.set_xticklabels(cpgPositions, rotation=90)
        ax.set_yticks([])
        if title:
            ax.set_title(title)
        radius = min(v_steps / 2.5, 0.05)
        h_spacings = self._buffer_spacings(h_spacings, radius * 2)

        if sort:
            cpgMatrix = self._sort_matrix_by_methylation(cpgMatrix)

        cpg_counter = 0
        for read, vspace in zip(cpgMatrix, v_spacings):
            ax.axhline(vspace, color="black", zorder=1)
            for cpg, hspace in zip(read, h_spacings):
                x = hspace
                y = vspace
                # only plot if a known value is provided
                if cpg == 1 or cpg==0:
                    if highlights is not None:
                        circle = plt.Circle((x, y), radius=radius, facecolor=self._get_color(cpg), 
                                            edgecolor=self._get_highlight_color(highlights[cpg_counter]), linewidth=3)
                        cpg_counter += 1
                    else:
                        circle = plt.Circle((x, y), radius=radius, facecolor=self._get_color(cpg), 
                                            edgecolor="black")
                    ax.add_artist(circle)
                else:
                    cpg_counter +=1
#         ax.axis("equal")
        return ax

    def plotCpGReadDataFrame(self, clustered_data_frame, title=None, figsize=(8, 8), sort=True):
        """
        Wrapper method to combine prepping a clustered data frame and plotting the matrix
        :param clustered_data_frame: #todo describe this
        :param figsize: tuple passed to matplotlib figsize (8,8)
        :return: matplotlib figure
        """
        prepped_data = self.prep_clustered_data_frame(clustered_data_frame)

        return self.plotCpGMatrix(prepped_data[0], prepped_data[1], title=title, figsize=figsize, sort=sort)