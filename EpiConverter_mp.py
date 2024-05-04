# from matchms.importing import load_from_mzxml
from matchms import set_matchms_logger_level
import csv
import os
import argparse
from datetime import datetime
from configparser import ConfigParser

from typing import Generator
import numpy as np
from pyteomics import mzxml
from matchms.importing.parsing_utils import parse_mzml_mzxml_metadata
from matchms.Spectrum import Spectrum
import timeit
import gc
import multiprocessing as mp
from tqdm import tqdm

from lxml.etree import XMLSyntaxError
import shutil

def load_from_mzxml(filename: str, ms_level: int = 2,
                    metadata_harmonization: bool = True) -> Generator[Spectrum, None, None]:
    """Load spectrum(s) from mzml file.

    This function will create ~matchms.Spectrum for every spectrum of desired
    ms_level found in a given MzXML file. For more extensive parsing options consider
    using the pyteomics package.

    Example:

    .. code-block:: python

        from matchms.importing import load_from_mzxml

        file_mzxml = "testdata.mzxml"
        spectrums = list(load_from_mzml(file_mzxml))

    Parameters
    ----------
    filename:
        Filename for mzXML file to import.
    ms_level:
        Specify which ms level to import. Default is 2.
    metadata_harmonization : bool, optional
        Set to False if metadata harmonization to default keys is not desired.
        The default is True.
    """
    for pyteomics_spectrum in mzxml.read(filename, dtype=dict):

        # print("pyteomics_spectrum", pyteomics_spectrum, type(pyteomics_spectrum))

        if ("ms level" in pyteomics_spectrum and pyteomics_spectrum["ms level"] == ms_level
                or "msLevel" in pyteomics_spectrum and pyteomics_spectrum["msLevel"] == ms_level):
            metadata = parse_mzml_mzxml_metadata(pyteomics_spectrum)

            # add additional metadata key-value pairs to metadata (added)
            if "precursorMz" in pyteomics_spectrum and "precursorScanNum" in pyteomics_spectrum["precursorMz"][0]:
                metadata["precursor_scannum"] = pyteomics_spectrum["precursorMz"][0]["precursorScanNum"]
            # else:
            #     metadata["precursor_scannum"] = None
            
            if "precursorMz" in pyteomics_spectrum and "activationMethod" in pyteomics_spectrum["precursorMz"][0]:
                metadata["activation_method"] = pyteomics_spectrum["precursorMz"][0]["activationMethod"]
            
            mz = np.asarray(pyteomics_spectrum["m/z array"], dtype="float")
            intensities = np.asarray(pyteomics_spectrum["intensity array"], dtype="float")

            if mz.shape[0] > 0:
                # Sort by mz (if not sorted already)
                if not np.all(mz[:-1] <= mz[1:]):
                    idx_sorted = np.argsort(mz)
                    mz = mz[idx_sorted]
                    intensities = intensities[idx_sorted]

                yield Spectrum(mz=mz, intensities=intensities, metadata=metadata,
                               metadata_harmonization=metadata_harmonization)

def creat_path(path):
    """Creat path+folder if not exist. Do nothing if path exists

    Parameters
    ----------
    path : str
        path + folder_name
    """
    isExists = os.path.exists(path)

    if not isExists:
        os.makedirs(path)
        # print("path generated:", path)
    else:
        # print("path exists:", path)
        pass

def find_all_dataset(input_path, dataset_format="mzXML"):
    """find all the dataset under input_path in the format of dataset_format

    Parameters
    ----------
    input_path : str
        input path of dataset
    dataset_format : str, optional
        format of dataset files, by default "mzXML"
    
    Returns
    -------
    dataset_list : list
        The list of dataset under input_path
    """
    dataset_list = []

    # find all the dataset under input_path
    dirs = os.listdir(input_path)
    # print("dirs", dirs)
    for dir in dirs:
        if dataset_format in dir:
            dataset_list.append(dir.split(".")[0])
    
    return dataset_list

def convert(dataset):
    """converting the mzMXL file to external files

    Parameters
    ----------
    dataset : str
        dataset name
    """
    set_matchms_logger_level("ERROR")

    fn_format = "mzXML"

    # read configuration from config_EpiConverter.ini
    cfg = ConfigParser()
    cfg.read("./config_EpiConverter.ini")
    cfg_dict = dict(cfg.items("config"))
    cfg_ms1_dict = dict(cfg.items("ms1"))
    cfg_ms2_dict = dict(cfg.items("ms2"))

    print("general configs:", cfg_dict)
    print("MS1 configs:", cfg_ms1_dict)
    print("MS2 configs:", cfg_ms2_dict)
    
    # general configs
    input_path = cfg_dict["input_path"]
    output_path = cfg_dict["output_path"]

    ion_mobility = False
    if "ion_mobility" in cfg_dict and cfg_dict["ion_mobility"] == "yes":
        ion_mobility = True

    # MS1 configs
    ion_injection_time_ms1 = cfg_ms1_dict["ion_injection_time"]
    instrument_type_ms1 = cfg_ms1_dict["instrument_type"]
    # MS2 configs
    ion_injection_time_ms2 = cfg_ms2_dict["ion_injection_time"]
    activation_type = cfg_ms2_dict["activation_type"]
    instrument_type_ms2 = cfg_ms2_dict["instrument_type"]
    acquisition_method = cfg_ms2_dict["acquisition_method"]


    # create the output path
    output_path_ms1 = output_path + "/MS1"
    output_path_ms2 = output_path + "/MS2"
    creat_path(output_path_ms1)
    creat_path(output_path_ms2)

    # check acquisition method
    if acquisition_method not in ["DIA", "DDA"]:
        print("Unrecognized acquisition method:", acquisition_method)
        exit()

    # starting processing each dataset
    start = timeit.default_timer()

    print("Processing " + dataset + "." + fn_format)

    file_name = input_path + "/" + dataset + "." + fn_format
    
    # check whether fake raw data exists. If it exists, skip processing
    if os.path.isfile(output_path + "/" + dataset + ".raw"):
        print("Fake raw data exists. Skip processing %s." % dataset)
        return
    
    # processing the MS1 file
    print("processing MS1...")

    if fn_format == "mzXML":
        spectrums_ms1_all_list = list(load_from_mzxml(file_name, ms_level=1))
    else:
        raise("Unsupported file format.")
    
    # get the current date and time
    now = datetime.now()
    dt_string = now.strftime("%m-%d-%Y %H:%M:%S")

    # print("spectrum:", spectrums_ms1_all_list[0])
    # print("spectrums_ms1_all_list[0].metadata:", spectrums_ms1_all_list[0].metadata)
    # print("spectrums_ms1_all_list[1].metadata:", spectrums_ms1_all_list[1].metadata)
    # print("spectrums_ms1_all_list[-1].metadata:", spectrums_ms1_all_list[-1].metadata)

    print("Generating the MS1 file...")
    with open(output_path_ms1 + "/" + dataset + ".MS1", "w") as f:
        # write header info to the ms1 result file
        f.write("H\tCreationDate\t" + dt_string + "\n")
        f.write("H\tExtractor\tEpiConverter\n")
        f.write("H\tExtractorVersion\t2.0\n")
        f.write("H\tDataType\tCentroid\n")

        # find the first scan_num
        scan_num_first = spectrums_ms1_all_list[0].metadata["scan_number"]
        # find the last scan_num
        scan_num_last = spectrums_ms1_all_list[-1].metadata["scan_number"]
        f.write("H\tFirstScan\t" + scan_num_first + "\n")
        f.write("H\tLastScan\t" + scan_num_last + "\n")

        # process each specturm in ms1
        for spectrum_ms1 in spectrums_ms1_all_list:
            f.write("S\t" + spectrum_ms1.metadata["scan_number"] + "\t" +  spectrum_ms1.metadata["scan_number"] + "\n")
            f.write("I\tRetTime\t" + "%.6f" % spectrum_ms1.metadata["retention_time"] + "\n")
            # f.write("I\tIonInjectionTime\t100\n")
            # f.write("I\tInstrumentType\tFTMS\n")
            f.write("I\tIonInjectionTime\t%s\n" % ion_injection_time_ms1)
            f.write("I\tInstrumentType\t%s\n" % instrument_type_ms1)

            peak_number = len(spectrum_ms1.peaks.mz)

            # if ino_mobility is no, directly write mz and corresponding intensity
            # if ino_mobility is yes, sum up intensities of similar/repetitive mz
            if not ion_mobility:
                for peak_idx in range(peak_number):
                    f.write("%.5f" % spectrum_ms1.peaks.mz[peak_idx] + " " + "%.1f" % spectrum_ms1.peaks.intensities[peak_idx] + "\n")
            else:
                if peak_number != 0:
                    cur_mz = spectrum_ms1.peaks.mz[0]
                    cur_intensity = 0
                
                for peak_idx in range(peak_number):
                    # find similar mz
                    if abs(cur_mz - spectrum_ms1.peaks.mz[peak_idx]) <= 0.001:
                        cur_intensity += spectrum_ms1.peaks.intensities[peak_idx]
                    # once find the next mz, write the current mz and intensity to the file and update the current mz and intensity
                    else:
                        f.write("%.5f" % cur_mz + " " + "%.1f" % cur_intensity + "\n")
                        cur_mz = spectrum_ms1.peaks.mz[peak_idx]
                        cur_intensity = spectrum_ms1.peaks.intensities[peak_idx]

                    # if it is the last peak, write the current mz and intensity to the file
                    if peak_idx == (peak_number - 1):
                        f.write("%.5f" % cur_mz + " " + "%.1f" % cur_intensity + "\n")

    # processing the MS2 file
    print("processing MS2...")
    if fn_format == "mzXML":
        spectrums_ms2_all_list = list(load_from_mzxml(file_name, ms_level=2))
    else:
        raise("Unsupported file format.")
    
    # print("spectrums_ms2_all_list[0].metadata", spectrums_ms2_all_list[0].metadata)
    # print("spectrums_ms2_all_list[-1].metadata", spectrums_ms2_all_list[-1].metadata)

    # Removing ms2 spectrums that exceed the scanNum range of ms1 for the new instrustment method
    if int(int(spectrums_ms1_all_list[0].metadata["scan_number"]) != 1):
        print("Removing ms2 spectrums that exceed the scanNum range of ms1")
        ms2_remove_idx_list = []

        ms1_scan_num_first = int(spectrums_ms1_all_list[0].metadata["scan_number"])
        ms1_scan_num_last = int(spectrums_ms1_all_list[-1].metadata["scan_number"])

        # for i, spectrum_ms2 in enumerate(spectrums_ms2_all_list):
        #     ms2_scan_num = int(spectrum_ms2.metadata["scan_number"])

        #     if ms2_scan_num < ms1_scan_num_first or ms2_scan_num > ms1_scan_num_last:
        #         spectrums_ms2_all_list.pop(i)

        # iterate over spectrums_ms1_all_list reversely
        for i in range(len(spectrums_ms2_all_list)-1, -1, -1):
            ms2_scan_num = int(spectrums_ms2_all_list[i].metadata["scan_number"])

            if ms2_scan_num < ms1_scan_num_first or ms2_scan_num > ms1_scan_num_last:
                spectrums_ms2_all_list.pop(i)
    else:
        print("Keep ms2 spectrums list")

    # print("spectrums_ms2_all_list[0].metadata", spectrums_ms2_all_list[0].metadata)
    # print("spectrums_ms2_all_list[-1].metadata", spectrums_ms2_all_list[-1].metadata)

    # add precursor_scannum to ms2 spectrums if it doesn't exist
    precursor_scannum_exist = True
    if "precursor_scannum" not in spectrums_ms2_all_list[0].metadata or spectrums_ms2_all_list[0].metadata["precursor_scannum"] is None:

        print("MS2 doesn't contain precursorScanNum. Calculating precursorScanNum...")
        precursor_scannum_exist = False

        ms1_scan_num_list = [int(spectrum_ms1.metadata["scan_number"]) for spectrum_ms1 in spectrums_ms1_all_list]
        ms1_scan_num_list.append(float("inf"))

        ms2_precursorScanNum_list = []

        # add precursorScanNum to each MS2 spectrum (time complexity O(n))
        if acquisition_method == "DIA":
            ms1_start_idx = 0
            ms1_end_idx = 1
            for spectrum_ms2 in spectrums_ms2_all_list:
                ms2_scan_num = int(spectrum_ms2.metadata["scan_number"])
                ms1_scan_num_start = ms1_scan_num_list[ms1_start_idx]
                ms1_scan_num_end = ms1_scan_num_list[ms1_end_idx]

                # print("Before processing", "ms2_scan_num:", ms2_scan_num, "ms1_scan_num_start:", ms1_scan_num_start, "ms1_scan_num_end:", ms1_scan_num_end)

                if ms2_scan_num > ms1_scan_num_start and ms2_scan_num < ms1_scan_num_end:
                    pass
                else:
                    ms1_start_idx += 1
                    ms1_end_idx += 1

                    ms1_scan_num_start = ms1_scan_num_list[ms1_start_idx]
                    ms1_scan_num_end = ms1_scan_num_list[ms1_end_idx]
                
                # print("After processing", "ms2_scan_num:", ms2_scan_num, "ms1_scan_num_start:", ms1_scan_num_start, "ms1_scan_num_end:", ms1_scan_num_end)

                assert(ms2_scan_num > ms1_scan_num_start and ms2_scan_num < ms1_scan_num_end)
                # spectrum_ms2.metadata.set("precursor_scannum", ms1_scan_num_start)
                # spectrum_ms2.metadata["precursor_scannum"] = ms1_scan_num_start
                # print("spectrum_ms2.metadata", spectrum_ms2.metadata)

                ms2_precursorScanNum_list.append(ms1_scan_num_start)
            
        elif acquisition_method == "DDA":
            ms1_start_idx = 0
            ms1_end_idx = 1
            for spectrum_ms2 in spectrums_ms2_all_list:
                ms2_scan_num = int(spectrum_ms2.metadata["scan_number"])
                ms1_scan_num_start = ms1_scan_num_list[ms1_start_idx]
                ms1_scan_num_end = ms1_scan_num_list[ms1_end_idx]

                # print("Before processing", "ms2_scan_num:", ms2_scan_num, "ms1_scan_num_start:", ms1_scan_num_start, "ms1_scan_num_end:", ms1_scan_num_end)

                # if ms2_scan_num > ms1_scan_num_start and ms2_scan_num < ms1_scan_num_end:
                #     pass
                # else:
                #     ms1_start_idx += 1
                #     ms1_end_idx += 1

                #     ms1_scan_num_start = ms1_scan_num_list[ms1_start_idx]
                #     ms1_scan_num_end = ms1_scan_num_list[ms1_end_idx]
                
                while ms2_scan_num <= ms1_scan_num_start or ms2_scan_num >= ms1_scan_num_end:
                    ms1_start_idx += 1
                    ms1_end_idx += 1

                    ms1_scan_num_start = ms1_scan_num_list[ms1_start_idx]
                    ms1_scan_num_end = ms1_scan_num_list[ms1_end_idx]
                
                # print("After processing", "ms2_scan_num:", ms2_scan_num, "ms1_scan_num_start:", ms1_scan_num_start, "ms1_scan_num_end:", ms1_scan_num_end)

                assert(ms2_scan_num > ms1_scan_num_start and ms2_scan_num < ms1_scan_num_end)
                # spectrum_ms2.metadata.set("precursor_scannum", ms1_scan_num_start)
                # spectrum_ms2.metadata["precursor_scannum"] = ms1_scan_num_start
                # print("spectrum_ms2.metadata", spectrum_ms2.metadata)

                ms2_precursorScanNum_list.append(ms1_scan_num_start)
        

    print("Generating the MS2 file...")
    with open(output_path_ms2 + "/" + dataset + ".ms2", "w") as f:
        # write header info to the ms2 result file
        f.write("H\tCreationDate\t" + dt_string + "\n")
        f.write("H\tExtractor\tEpiConverter\n")
        f.write("H\tExtractorVersion\t2.0\n")
        f.write("H\tComments\tCreated by Chenfeng, 2023\n")
        f.write("H\tDataType\tCentroid\n")
        
        for i, spectrum_ms2 in enumerate(spectrums_ms2_all_list):
            f.write("S\t" + spectrum_ms2.metadata["scan_number"] + "\t" + spectrum_ms2.metadata["scan_number"] + "\t" + str(int(spectrum_ms2.metadata["precursor_mz"])) + "\n")
            f.write("I\tNumberOfPeaks\t" + str(len(spectrum_ms2.peaks.mz)) +"\n")
            f.write("I\tRetTime\t%.6f\n" % spectrum_ms2.metadata["retention_time"])
            # f.write("I\tIonInjectionTime\t10\n")
            f.write("I\tIonInjectionTime\t%s\n" % ion_injection_time_ms2)

            # if "activation_method" in spectrum_ms2.metadata:
            #     f.write("I\tActivationType\t" + spectrum_ms2.metadata["activation_method"] + "\n")
            # else:
            #     f.write("I\tActivationType\t" + activation_type + "\n")
            
            f.write("I\tActivationType\t" + activation_type + "\n")
            # f.write("I\tInstrumentType\tFTMS\n")
            f.write("I\tInstrumentType\t%s\n" % instrument_type_ms2)

            # precursorScanNum
            if not precursor_scannum_exist:
                f.write("I\tPrecursorScan\t" + str(ms2_precursorScanNum_list[i]) + "\n")
            else:
                f.write("I\tPrecursorScan\t" + spectrum_ms2.metadata["precursor_scannum"] + "\n")
            
            if acquisition_method == "DIA":
                f.write("I\tActivationCenter\t" + str(spectrum_ms2.metadata["precursor_mz"]) + "\n")
                f.write("I\tMonoiosotopicMz\t0\n")
                # f.write("Z\t??")
            elif acquisition_method == "DDA":
                f.write("I\tActivationCenter\t%.2f\n" % spectrum_ms2.metadata["precursor_mz"])
                f.write("I\tMonoiosotopicMz\t%f\n" % spectrum_ms2.metadata["precursor_mz"])

            peak_number = len(spectrum_ms2.peaks.mz)

            # if ino_mobility is no, directly write mz and corresponding intensity
            # if ino_mobility is yes, sum up intensities of similar/repetitive mz
            if not ion_mobility:
                for peak_idx in range(peak_number):
                    f.write("%.5f %.1f\n" % (spectrum_ms2.peaks.mz[peak_idx], spectrum_ms2.peaks.intensities[peak_idx]))
            else:
                if peak_number != 0:
                    cur_mz = spectrum_ms2.peaks.mz[0]
                    cur_intensity = 0
                
                for peak_idx in range(peak_number):
                    # find similar mz
                    if abs(cur_mz - spectrum_ms2.peaks.mz[peak_idx]) <= 0.001:
                        cur_intensity += spectrum_ms2.peaks.intensities[peak_idx]
                    # once find the next mz, write the current mz and intensity to the file and update the current mz and intensity
                    else:
                        f.write("%.5f" % cur_mz + " " + "%.1f" % cur_intensity + "\n")
                        cur_mz = spectrum_ms2.peaks.mz[peak_idx]
                        cur_intensity = spectrum_ms2.peaks.intensities[peak_idx]

                    # if it is the last peak, write the current mz and intensity to the file
                    if peak_idx == (peak_number - 1):
                        f.write("%.5f" % cur_mz + " " + "%.1f" % cur_intensity + "\n")

    # delete ms data
    del spectrums_ms1_all_list
    del spectrums_ms2_all_list
    gc.collect()

    # create a fake file at output_path
    # os.system("touch " + output_path + "/" + dataset + ".raw")
    with open(output_path + "/" + dataset + ".raw", "w") as fp:
        pass

    end = timeit.default_timer()
    print("Running time: %s Seconds"%(end-start))


def filter_datafile(dataset):
    """This function is to filter data files with syntax errors (lxml.etree.XMLSyntaxError) by finding and moving them to a other dir

    Parameters
    ----------
    dataset : str
        dataset name
    """

    set_matchms_logger_level("ERROR")

    fn_format = "mzXML"

    # read configuration from config_EpiConverter.ini
    cfg = ConfigParser()
    cfg.read("./config_EpiConverter.ini")
    cfg_dict = dict(cfg.items("config"))
    cfg_ms1_dict = dict(cfg.items("ms1"))
    cfg_ms2_dict = dict(cfg.items("ms2"))

    # general configs
    input_path = cfg_dict["input_path"]
    output_path = cfg_dict["output_path"]

    # create path to place data files with syntax errors
    creat_path(input_path + "/" + "wrong_data")

    file_name = input_path + "/" + dataset + "." + fn_format

    # check whether fake raw data exists. If it exists, skip filtering
    if os.path.isfile(output_path + "/" + dataset + ".raw"):
        print("Fake raw data exists. Skip filtering %s." % dataset)
        return

    # check whether the datafile could be recognized
    try:
        if fn_format == "mzXML":
            spectrums_ms1_all_list = list(load_from_mzxml(file_name, ms_level=1))
            spectrums_ms2_all_list = list(load_from_mzxml(file_name, ms_level=2))

            # delete ms data
            del spectrums_ms1_all_list
            del spectrums_ms2_all_list
            gc.collect()
        else:
            raise("Unsupported file format.")
    except XMLSyntaxError:
        new_file_name = input_path + "/" + "wrong_data" + "/" + dataset + "." + fn_format
        print("Syntax errors are found in " + dataset + ", moving it to " + new_file_name)
        shutil.move(file_name, new_file_name)



if __name__ == "__main__":
    # set_matchms_logger_level("ERROR")

    # input_path = "./data/EpiConverter"
    # output_path = "./result/EpiConverter"

    fn_format = "mzXML"

    # read configuration from config_EpiConverter.ini
    cfg = ConfigParser()
    cfg.read("./config_EpiConverter.ini")
    cfg_dict = dict(cfg.items("config"))
    cfg_ms1_dict = dict(cfg.items("ms1"))
    cfg_ms2_dict = dict(cfg.items("ms2"))

    # general configs
    input_path = cfg_dict["input_path"]
    core_num = int(cfg_dict["core_num"])
    if "filter_dataset" in cfg_dict.keys():
        filter_dataset = eval(cfg_dict["filter_dataset"].lower().capitalize())
    else:
        filter_dataset = False

    # find all the dataset under input_path and filter the dataset with syntax errors
    if filter_dataset:
        print("dataset filter is enabled.")
        dataset_list = find_all_dataset(input_path, fn_format)
        for dataset in tqdm(dataset_list, desc="Filtering Dataset"):
            # finding and moving dataset with syntax errors to another dir
            filter_datafile(dataset)
    else:
        print("dataset filter is disabled.")


    # find all the dataset under input_path again
    dataset_list = find_all_dataset(input_path, fn_format)
    print("dataset_list:", dataset_list)
    
    print("Start parallel computing with %d cores" % core_num)
    
    # for dataset in dataset_list:
    #     convert(dataset)

    # parallel converting
    # pool = mp.Pool(processes=core_num)
    # pool.map(convert, dataset_list)

    with mp.Pool(processes=core_num) as pool:
        for _ in tqdm(pool.imap_unordered(convert, dataset_list), desc="Dataset Processing", total=len(dataset_list)):
            pass


    print("Converting Finished!")
    










        

