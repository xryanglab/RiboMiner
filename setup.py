#!/usr/bin/env python
# -*- coding:UTF-8 -*-
import os, sys
from setuptools import setup
from RiboMiner import __version__

if sys.version_info.major != 3:
    sys.exit("RiboAnalyzer can only be used with Python 3. You are currently "
             "running Python %d." % sys.version_info.major)

with open("README.md", "r",encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='RiboMiner',
    version = __version__,
    description = 'A intergrated tool for downstream analysis of ribosome profiling data.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="Analysis for ribosome profiling data.",
    url="https://github.com/sherkinglee/RiboAnalyzer",
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
            'cAI=RiboAnalyzer.cAI:main',
            'EnrichmentAnalysis=RiboAnalyzer.EnrichmentAnalysis:main',
            'EnrichmentAnalysisForSingleTrans=RiboAnalyzer.EnrichmentAnalysisForSingleTrans:main',
            'enrichmentMeanDensity=RiboAnalyzer.enrichmentMeanDensity:main',
            'ExtractSequenceCenteredOnAPosition=RiboAnalyzer.ExtractSequenceCenteredOnAPosition:main',
            'GCContent=RiboAnalyzer.GCContent:main',
            'GetProteinCodingSequence=RiboAnalyzer.GetProteinCodingSequence:main',
            'GetUTRSequences=RiboAnalyzer.GetUTRSequences:main',
            'hydropathyCharge=RiboAnalyzer.hydropathyCharge:main',
            'LengthDistribution=RiboAnalyzer.LengthDistribution:main',
            'MergeSampleDensitys=RiboAnalyzer.MergeSampleDensitys:main',
            'MetageneAnalysis=RiboAnalyzer.MetageneAnalysis:main',
            'MetageneAnalysisForTheWholeRegions=RiboAnalyzer.MetageneAnalysisForTheWholeRegions:main',
            'OutputTranscriptInfo=RiboAnalyzer.OutputTranscriptInfo:main',
            'PausingScore=RiboAnalyzer.PausingScore:main',
            'Periodicity=RiboAnalyzer.Periodicity:main',
            'PlotEnrichmentRatio=RiboAnalyzer.PlotEnrichmentRatio:main',
            'PlotGCContent=RiboAnalyzer.PlotGCContent:main',
            'PlotHydropathyCharge=RiboAnalyzer.PlotHydropathyCharge:main',
            'PlotMetageneAnalysis=RiboAnalyzer.PlotMetageneAnalysis:main',
            'PlotMetageneAnalysisForTheWholeRegions=RiboAnalyzer.PlotMetageneAnalysisForTheWholeRegions:main',
            'PlotPolarity=RiboAnalyzer.PlotPolarity:main',
            'PlotRiboDensityAroundTriAAMotifs=RiboAnalyzer.PlotRiboDensityAroundTriAAMotifs:main',
            'PlotRiboDensityAtEachKindAAOrCodon=RiboAnalyzer.PlotRiboDensityAtEachKindAAOrCodon:main',
            'PolarityCalculation=RiboAnalyzer.PolarityCalculation:main',
            'ProcessPausingScore=RiboAnalyzer.ProcessPausingScore:main',
            'RiboDensityAroundTripleteAAMotifs=RiboAnalyzer.RiboDensityAroundTripleteAAMotifs:main',
            'RiboDensityAtEachKindAAOrCodon=RiboAnalyzer.RiboDensityAtEachKindAAOrCodon:main',
            'RiboDensityAtEachPosition=RiboAnalyzer.RiboDensityAtEachPosition:main',
            'RiboDensityForSpecificRegion=RiboAnalyzer.RiboDensityForSpecificRegion:main',
            'RiboDensityOfDiffFrames=RiboAnalyzer.RiboDensityOfDiffFrames:main',
            'RPFdist=RiboAnalyzer.RPFdist:main',
            'StatisticReadsOnDNAsContam=RiboAnalyzer.StatisticReadsOnDNAsContam:main',
            'tAI=RiboAnalyzer.tAI:main',
            'tAIPlot=RiboAnalyzer.tAIPlot:main',
            'cAIPlot=RiboAnalyzer.cAIPlot:main',
        ],

    },
        )