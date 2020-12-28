#!/usr/bin/env python
# -*- coding:UTF-8 -*-
import os, sys
from setuptools import setup
from RiboMiner import __version__

if sys.version_info.major != 3:
    sys.exit("RiboMiner can only be used with Python 3. You are currently "
             "running Python %d." % sys.version_info.major)

with open("README.md", "r",encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='RiboMiner',
    version = __version__,
    description = 'A python toolset for mining multi-dimensional features of the translatome with ribosome profiling data',
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="Analysis for ribosome profiling data.",
    url="https://github.com/xryanglab/RiboMiner",
    author = 'Li Fajin',
    author_email = 'sherkinglee@gmail.com',
    license='GPLv3.0',
    packages=['RiboMiner','data'],
    install_requires=[
                      'matplotlib>=2.1.0',
                      'numpy>=1.16.4',
                      'pandas>=0.24.2',
                      'pysam>=0.15.2',
                      'scipy>=1.1.0',
                      'seaborn>=0.8.1',
                      'biopython>=1.70',
                      'scipy>=1.1.0',
                      'HTSeq',
                      'RiboCode>=1.2.10',
                      'pysamstats',
                      ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Environment :: Console',
        'Programming Language :: Python :: 3',
        'Operating System :: Unix',
        ],
    python_requires='>=3.6',
    entry_points={

        'console_scripts': [
            #name determined the name of cmd line direct call
            'cAI=RiboMiner.cAI:main',
            'EnrichmentAnalysis=RiboMiner.EnrichmentAnalysis:main',
            'EnrichmentAnalysisForSingleTrans=RiboMiner.EnrichmentAnalysisForSingleTrans:main',
            'enrichmentMeanDensity=RiboMiner.enrichmentMeanDensity:main',
            'ExtractSequenceCenteredOnAPosition=RiboMiner.ExtractSequenceCenteredOnAPosition:main',
            'GCContent=RiboMiner.GCContent:main',
            'GetProteinCodingSequence=RiboMiner.GetProteinCodingSequence:main',
            'GetUTRSequences=RiboMiner.GetUTRSequences:main',
            'hydropathyCharge=RiboMiner.hydropathyCharge:main',
            'LengthDistribution=RiboMiner.LengthDistribution:main',
            'MergeSampleDensitys=RiboMiner.MergeSampleDensitys:main',
            'MetageneAnalysis=RiboMiner.MetageneAnalysis:main',
            'MetageneAnalysisForTheWholeRegions=RiboMiner.MetageneAnalysisForTheWholeRegions:main',
            'OutputTranscriptInfo=RiboMiner.OutputTranscriptInfo:main',
            'PausingScore=RiboMiner.PausingScore:main',
            'Periodicity=RiboMiner.Periodicity:main',
            'PlotEnrichmentRatio=RiboMiner.PlotEnrichmentRatio:main',
            'PlotGCContent=RiboMiner.PlotGCContent:main',
            'PlotHydropathyCharge=RiboMiner.PlotHydropathyCharge:main',
            'PlotMetageneAnalysis=RiboMiner.PlotMetageneAnalysis:main',
            'PlotMetageneAnalysisForTheWholeRegions=RiboMiner.PlotMetageneAnalysisForTheWholeRegions:main',
            'PlotPolarity=RiboMiner.PlotPolarity:main',
            'PlotRiboDensityAroundTriAAMotifs=RiboMiner.PlotRiboDensityAroundTriAAMotifs:main',
            'PlotRiboDensityAtEachKindAAOrCodon=RiboMiner.PlotRiboDensityAtEachKindAAOrCodon:main',
            'PolarityCalculation=RiboMiner.PolarityCalculation:main',
            'ProcessPausingScore=RiboMiner.ProcessPausingScore:main',
            'RiboDensityAroundTripleteAAMotifs=RiboMiner.RiboDensityAroundTripleteAAMotifs:main',
            'RiboDensityAtEachKindAAOrCodon=RiboMiner.RiboDensityAtEachKindAAOrCodon:main',
            'RiboDensityAtEachPosition=RiboMiner.RiboDensityAtEachPosition:main',
            'RiboDensityForSpecificRegion=RiboMiner.RiboDensityForSpecificRegion:main',
            'RiboDensityOfDiffFrames=RiboMiner.RiboDensityOfDiffFrames:main',
            'RPFdist=RiboMiner.RPFdist:main',
            'StatisticReadsOnDNAsContam=RiboMiner.StatisticReadsOnDNAsContam:main',
            'tAI=RiboMiner.tAI:main',
            'tAIPlot=RiboMiner.tAIPlot:main',
            'cAIPlot=RiboMiner.cAIPlot:main',
            'ModifyHTseq=RiboMiner.ModifyHTseq:main',
            'ReadsLengthOfSpecificRegions=RiboMiner.ReadsLengthOfSpecificRegions:main',
            'CoverageOfEachTrans=RiboMiner.CoverageOfEachTrans:main',
            'PlotTransCoverage=RiboMiner.PlotTransCoverage:main',
        ],

    },
        )