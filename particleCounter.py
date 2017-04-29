import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from skimage.segmentation import mark_boundaries
from skimage.filters import threshold_otsu
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
import os
from PIL import Image

SEGMENTFILENAME = 'segmented.pdf'
FILEPATH = os.path.abspath('example2.tif')
SEM_SCALE = 0.005  # um per pixel
SEM_SCALE = 0.0227794
PIXEL_AREA = SEM_SCALE ** 2
HISTOGRAM_FILENAME = 'histogram.pdf'
histogram_path = 'size_distribution'
histogram_title = 'Histogram'
histogram_width = 4
histogram_height = 4
histogram_binwidth = 0.025
histogram_xlimmax = 0.5
MULTIPLIER = 3.2
MAX_DIAMETER = 0.6     #um
MIN_DIAMETER = 0.03    #um
BLURRING = 30
imageAuto = True


def treat_image_auto(image_path, pixel_size):
    imagea = Image.open(image_path)

    original_gray = np.asarray(imagea.convert('L'), dtype=float)
    blurred = ndi.gaussian_filter(original_gray, BLURRING)
    difference = original_gray - (MULTIPLIER * blurred)
    threshold = original_gray >= threshold_otsu(original_gray)
    imagetolabel = ((difference > threshold) * 255).astype('uint8')

    local_maxi = peak_local_max(imagetolabel, indices=False, footprint=np.ones([21, 21]), labels=threshold)
    labels, count = ndi.label(local_maxi)

    remove_index = []
    particle_diameters = []
    for particle_index in range(count):
        num_pixels = (labels == particle_index).sum()
        print 'num_pixels: ', num_pixels
        area = num_pixels * pixel_size
        p_diameter = 2 * np.sqrt(area / np.pi)
        print 'diameter: ', p_diameter
        if p_diameter > MAX_DIAMETER:
            print 'not counting (high size): ', particle_index
            remove_index.append(particle_index)
        elif p_diameter < MIN_DIAMETER:
            print 'not counting (low size): ', particle_index
        else:
            particle_diameters.append(p_diameter)

    print 'previous labels: ', len(labels)
    labels = np.delete(labels, remove_index, 0)
    print 'final labels: ', len(labels)

    print 'previous number Ps: ', count
    count = len(particle_diameters)
    print 'final number Ps: ', count

    return count, labels, particle_diameters


def treat_image_ImageJ(image_path, pixel_size):
    """Function to count particles from images treated like this in ImageJ:
        1. Make binary
        2. Gaussian blur
        3. Invert
    """
    imagea = Image.open(image_path)

    original_gray = (np.asarray(imagea, dtype=float))
    blurred = ndi.gaussian_filter(original_gray, BLURRING)
    difference = original_gray - (MULTIPLIER * blurred)
    threshold = original_gray >= threshold_otsu(original_gray)
    imagetolabel = ((difference > threshold) * 255).astype('uint8')

    local_maxi = peak_local_max(imagetolabel, indices=False,
                                 footprint=np.ones([21, 21]),
                                 labels=threshold)
    labels, count = ndi.label(local_maxi)

    remove_index = []
    particle_diameters = []
    for particle_index in range(count):
        num_pixels = (labels == particle_index).sum()
        print 'num_pixels: ', num_pixels
        area = num_pixels * pixel_size
        p_diameter = 2 * np.sqrt(area / np.pi)
        print 'diameter: ', p_diameter
        if p_diameter > MAX_DIAMETER:
            print 'not counting (high size): ', particle_index
            remove_index.append(particle_index)
        elif p_diameter < MIN_DIAMETER:
            print 'not counting (low size): ', particle_index
        else:
            particle_diameters.append(p_diameter)

    print 'previous labels: ', len(labels)
    labels = np.delete(labels, remove_index, 0)
    print 'final labels: ', len(labels)

    print 'previous number Ps: ', count
    count = len(particle_diameters)
    print 'final number Ps: ', count

    return count, labels, particle_diameters


def generate_particle_segments_graphics (labels):
    pp = PdfPages(SEGMENTFILENAME)
    fig = plt.figure()

    #plt.imshow(labels)
    #plt.show()

    super_marked = 1 - mark_boundaries(labels != 0, labels, color=(1,0,0), outline_color=(1,0,0))
    plt.imshow(super_marked[::-1,:])
    plt.xlim((0, labels.shape[1]))
    plt.ylim((0, labels.shape[0]))

    plt.tight_layout()

    #plt.show()

    pp.savefig()
    pp.close()

def create_histogram_plot(histogram_particle_diameters):
    pp = PdfPages(HISTOGRAM_FILENAME)

    fig = plt.figure()

    fig.set_size_inches((histogram_width, histogram_height))

    bin_size = histogram_binwidth; min_edge = 0; max_edge = max(histogram_particle_diameters)
    #bin_size = 0.05; min_edge = 0; max_edge = HISTOGRAM_MAX_WIDTH
    N = (max_edge-min_edge)/bin_size; Nplus1 = N + 1
    bin_list = np.linspace(min_edge, max_edge, Nplus1)

    n, bins = np.histogram(histogram_particle_diameters, bins=bin_list)

    probability_mass_function = 100 * n / sum(n)
    plt.bar(bins[0:-1], probability_mass_function, width=bin_size)

    plt.title(histogram_title)
    plt.xlim(0, histogram_xlimmax)
    plt.xlabel(r"Diameter [$\mu$m]")
    plt.ylabel("Frequency (%)")

    plt.tight_layout()
    #plt.savefig(HISTOGRAM_FILENAME)

    pp.savefig()
    pp.close()


#Create textfile with diameters list
def create_textfile_diameters(particle_diameters):
    textfile = open('diameters.txt', 'w')
    textfile.write("diameters: %s" % particle_diameters)
    textfile.close()


# Calculate median size of the particles
def calculate_median(histogram_particle_diameters):
    median_size = np.median(histogram_particle_diameters)
    mean_size = np.mean(histogram_particle_diameters)
    #NOTE: Area calculated for spherical particles
    mean_area = 3.141592 * (mean_size / 2) ** 2
    print '{0} um median size'.format(median_size)
    print '{0} um mean size'.format(mean_size)
    print '{0} um mean area'.format(mean_area)


if __name__ == '__main__':
    histogram_particle_diameters = []
    if imageAuto:
        num_particles, labeled_image, particle_diameters = treat_image_auto(FILEPATH, PIXEL_AREA)
    else:
        num_particles, labeled_image, particle_diameters = treat_image_ImageJ(FILEPATH, PIXEL_AREA)

    histogram_particle_diameters.extend(particle_diameters)

    #Create segmented figure
    generate_particle_segments_graphics(labeled_image)

    #Create textfile with diameters list
    create_textfile_diameters(particle_diameters)

    #Calculate median size of the particles
    calculate_median(histogram_particle_diameters)

    #Call function to create histogram figure
    create_histogram_plot(histogram_particle_diameters)

    #TODO: calculate density in global image
    #TODO: calculate density in specific areas
    #TODO: make an approximated profile surface image (3D) from NPs (substrate->flat)

