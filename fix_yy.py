#!/bin/env python3
#-*- coding: utf-8 -*-

"""
Mission: rename files:

+ 1. Sample sheet

  - to: YY168_20220906_DLJ_CnT
  - from:

+ 2. to company

  - to: YY194_20221103_生物物理所俞洋组_派森诺
  - from: 【派森诺样品信息单】文库类样品20221103_生物物理所_俞洋组_YY194—寄送地：上海

+ 3. QC report:

  - to: 文库质检报告_SP220908575_派森诺.pdf
  - from: 派森诺生物——SP220908575 王明老师文库质检报告.pdf

  - to: 文库质检报告_SP220908575_派森诺.pdf
  - from: SP220908575.pdf
"""


import os
import sys
import re
import fnmatch
import shutil
import warnings
import time
import pandas as pd


def update_obj(obj, d, force=True, remove=False):
    """
    Update the object, by dict
    Parameters:
    -----------
    obj : object
        an object
    d : dict
        a dict, save key:value pairs for object attribute:value
    force : bool
        update exists attributes to object, default: True
    remove : bool
        remove exists attributes from object, default: False
    """
    if remove is True:
        for k in obj.__dict__:
            delattr(obj, k)
    # add attributes
    if isinstance(d, dict):
        for k, v in d.items():
            if not hasattr(obj, k) or force:
                setattr(obj, k, v)
    return obj


def list_dir(x, recursive=False, include_dirs=False):
    out = []
    if isinstance(x, str):
        if os.path.isdir(x):
            n = 0 # levels
            for (root, d, f) in os.walk(x):
                dirs = [os.path.join(root, i) for i in d]
                files = [os.path.join(root, i) for i in f]
                out += files
                if include_dirs:
                    out += dirs
                if not recursive:
                    break # first level
        else:
            print(f'list_dir() skipped, x not a directory: {x}')
    else:
        print('list_dir() skipped, x expect str, got {type(x).__name__}')
    return sorted(out)


def list_file(path='.', pattern='*', recursive=False,
    include_dirs=False):
    files = list_dir(path, recursive, include_dirs)
    files = [f for f in files if fnmatch.fnmatch(os.path.basename(f), pattern)]
    return sorted(files)


class YYdir(object):
    """
    rename files in YY dir
    1. sample_sheet:
    2. to_company:
    3. gel_report:
    4. qc_report_summary:
    5. qc_report:
    """
    def __init__(self, yy_dir, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.yy_dir = yy_dir.strip('/')
        self.init_args()


    def init_args(self):
        self.lab = '生物物理所俞洋组'
        self.co = ['派森诺', '吉因加', '诺禾'] # companies
        self.yy = self.get_yy()
        self.date = self.get_date()
        self.xlsx_list = list_file(self.yy_dir, "*.xlsx")
        self.pptx_list = list_file(self.yy_dir, "*.pptx", recursive=True)
        self.pdf_list = list_file(self.yy_dir, "*.pdf", recursive=True)
        self.is_valid_yy_dir = len(self.xlsx_list) > 0
        # self.df_ss = self.load_ss(ss)
        if self.is_valid_yy_dir:
            self.files = self.sanitize_xlsx()
            self.files.update(self.sanitize_pptx())
            self.files.update(self.sanitize_pdf())
        else:
            self.files = {}


    def get_yy(self):
        p = re.compile('^((Y\w)(\d{2,3}))')
        m = p.search(os.path.basename(self.yy_dir))
        if isinstance(m, re.Match): # fix width
            yy = f'{m.group(2)}{int(m.group(3)):03d}'
        else:
            yy = os.path.basename(self.yy_dir).split('_')[0] #
            print(f'unknown directory: {self.yy_dir}')
        return yy


    def get_date(self):
        """
        Extract date from filename, 8-digits, startswith('20')
        """
        p = re.compile('(20\d{6})')
        m = p.search(os.path.basename(self.yy_dir))
        if isinstance(m, re.Match):
            dt = m.group(1)
        else:
            dt = '20220101'
        return dt


    def is_ss(self, x):
        """
        is sample_sheet file:
        1. YY168_20220906_DLJ_CnT.xlsx
        2. load_ss
        """
        if isinstance(x, str):
            p1 = re.compile('(派森诺|吉因加|诺禾|信息单|文库|生物物理|俞洋)') # not ss
            m1 = p1.search(os.path.basename(x)) # not ss
            p2 = re.compile('^Y[A-Z](\d{2,3}).*.xlsx')
            m2 = p2.search(os.path.basename(x))
            if isinstance(m1, re.Match):
                out = False
            elif isinstance(m2, re.Match):
                out = True
            else:
                df = self.load_ss(x)
                out = len(df) > 0
        else:
            print(f'unknown sample_sheet: {x}')
            out = False
        return out


    def sanitize_pdf(self):
        """
        POD: 派森诺生物——SP221105386 王明老师文库质检报告.pdf, SP221105386.pdf
        Gene+: 文库质控报告20221004-生物物理所YY组自建库测序技术服务.pdf
        rename:
        YY168_20220906_文库质检报告_01.pdf
        """
        d = {}
        for i, f in enumerate(self.pdf_list):
            p1 = re.compile('SP\d+|LPL\d+|派森诺') # POD
            p2 = re.compile('YY组自建库|MGI|吉因加') # Gene+
            p3 = re.compile('Novogene|诺禾', flags=re.IGNORECASE)
            p4 = re.compile('华大', flags=re.IGNORECASE)
            # f = os.path.basename(f)
            if p1.search(f):
                cc = '派森诺'
                m2 = re.search('(SP\d+|LPL\d+)', f)
                if isinstance(m2, re.Match):
                    cc = f'{m2.group(1)}_{cc}'
            elif p2.search(f):
                cc = '吉因加'
            elif p3.search(f):
                cc = '诺禾致源'
            elif p4.search(f):
                cc = '华大基因'
            else:
                cc = '测序公司'
            fx = f'{self.yy}_{self.date}_质检报告_{cc}_{i+1:02d}.pdf'
            d.update({f:fx})
        return d


    def sanitize_pptx(self):
        """
        gel report
        """
        d = {}
        if len(self.pptx_list) == 1:
            fx = f'{self.yy}_胶图.pptx'
            # f = os.path.basename(self.pptx_list[0])
            f = self.pptx_list[0]
            d.update({f:fx})
        elif len(self.pptx_list) > 1:
            for i, f in enumerate(self.pptx_list):
                # f = os.path.basename(f)
                fx = f'{self.yy}_胶图_{i+1:02d}.pptx'
                d.update({f:fx})
        else:
            pass
        return d


    def sanitize_xlsx(self):
        """
        rename excel files
        sample_sheet: YY168_20220906_DLJ_CnT.xlsx
        to_company: YY168_20220906_生物物理所俞洋组_派森诺.xlsx
        """
        d = {}
        # 1. xlsx sample_sheet, sample_to_company
        ss = [i for i in self.xlsx_list if self.is_ss(i)]
        sx = [i for i in self.xlsx_list if i not in ss]
        if len(ss) == 0:
            self.is_valid_yy_dir = False
            return None
        fx = self.fix_ss_name(ss[0])
        if len(ss) == 1:
            d.update({ss[0]:fx})
        else:
            for i,f in enumerate(ss):
                fx = self.fix_ss_name(f)
                fx = fx.replace('.xlsx', f'_{i+1:02d}.xlsx')
                d.update({f:fx})
        # rename sx: to company
        p = re.compile('(派森诺|吉因加|诺禾|文库信息单|华大基因)') # POD
        if len(sx) > 1:
            pass
        elif len(sx) == 1:
            # f = os.path.basename(sx[0])
            # f = sx[0]
            m = p.search(sx[0])
            if isinstance(m, re.Match):
                co = m.group(1)
                co = '诺禾' if co == '文库信息单' else co
            else:
                co = '测序公司'
            fx = f'{self.yy}_{self.date}_{self.lab}_{co}.xlsx'
            d.update({sx[0]:fx})
        else:
            for i,f in sx:
                # p = re.compile('(派森诺|吉因加)') # POD
                m = p.search(f)
                if isinstance(m, re.Match):
                    co = m.group(1)
                    co = '诺禾' if co == '文库信息单' else co
                else:
                    co = '测序公司'
                fx = f'{self.yy}_{self.date}_{self.lab}_{co}_{i+1:02d}.xlsx'
                d.update({f:fx})
        return d


    def load_ss(self, x):
        """
        sample sheet:
        """
        c = ['Lib_number*', 'Lib_user', 'Lib_type*']
        try:
            warnings.simplefilter(action='ignore', category=UserWarning)
            df = pd.read_excel(x, sheet_name='sample_sheet')
            df = df.dropna(thresh=2)
            df = df[c]
            df.columns = ['yy', 'user', 'lib_type']
            # df['name'] = self.sanitize(df['name'].to_list()) # sanitize
        except:
            df = pd.DataFrame(columns=['yy', 'user', 'lib_type'])
        return df


    def fix_ss_name(self, x):
        """
        input: YY168_20220906_DLJ_CnT.xlsx
        fix:
        ## lib_type
        mRNAseq  : RNA
        smRNAseq : smRNA
        ChIPseq  : ChIP
        ATAC-seq : ATAC
        ATACseq  : ATAC
        RiboSeq  : Ribo
        Others   : ...
        ## user
        JG : GJQ
        LW : LWW
        """
        try:
            df = self.load_ss(x)
            yy = df.loc[0, 'yy']
            user = df.loc[0,'user']
            user = re.sub('[0-9]+', '', user)
            user = user.strip()
            lib_type = df.loc[0, 'lib_type']
            ## fix user
            du = {'JG': 'GJQ', 'LW': 'LWW'}
            user = du.get(user, user)
            ## fix lib_type
            # p1 = re.compile('(YY\d{2,3})_') # yy
            # m1 = p1.search(x) #
            p2 = re.compile('([A-Za-z]+).xlsx') # lib_type
            m2 = p2.search(os.path.basename(x)) #
            if isinstance(m2, re.Match):
                m2g = m2.group(1)
                tt = [
                    'ATAC', 'ATACseq', 'ChIP', 'ChIPSeq', 'ChrRNA', 'CLIP', 'CLIPseq', 
                    'CnR', 'CnT', 'CUTAC', 'nanobody', 'MiniATAC', 'PROseq', 'Riboseq', 
                    'RNAseq', 'scifi', 'sgRNA', 'SHERRY', 'Smart', 'smRNA', 'TTseq',
                    ]
                if m2g in tt or lib_type == 'Others':
                    lib_type = m2g
            # if lib_type == 'Others': 
            #     p2 = re.compile('-|_([A-Za-z]).xlsx') # lib_type
            #     m2 = p2.search(os.path.basename(x)) #
            #     lib_type = m2.group(1)
            p3 = re.compile('[^\w]|(seq)', flags=re.IGNORECASE)
            lib_type = p3.sub('', lib_type)
            return f'{self.yy}_{self.date}_{user}_{lib_type}.xlsx'
        except:
            print('!AAAA-1', x)


def fix_yy(x, out_dir):
    """
    fix files in YY dir
    """
    if isinstance(x, str):
        x = x.strip('/')
        yy = YYdir(x)
        if yy.is_valid_yy_dir:
            dest_dir = os.path.join(out_dir, os.path.basename(x))
            if not os.path.exists(dest_dir):
                os.makedirs(dest_dir)
            for k,v in yy.files.items():
                dest_file = os.path.join(dest_dir, v)
                try:
                    if not os.path.exists(dest_file):
                        shutil.copy(k, dest_file)
                        # shutil.copy(k, dest_file)
                except:
                    print('could not copy file: {}'.format(dest_file))
    else:
        print(f'x expect str, got {type(x).__name__}')


def fix_yy2(src, dest, n_max=0, overwrite=False):
    # level-1, src is yy_dir
    yy = YYdir(src)
    if yy.is_valid_yy_dir:
        fix_yy(src, dest)
    else:
        d_list = list_file(src, 'Y*', include_dirs=True, recursive=False)
        if n_max > 0:
            d_list = d_list[-n_max:]
        if len(d_list) > 0:
            for i,f in enumerate(d_list):
                print(f'[{i+1:03d}/{len(d_list):03d}] - {os.path.basename(f)} {" ":>40s}', end='\r')
                # time.sleep(1)
                fix_yy(f, dest)
                # dest_yy = os.path.join(dest, os.path.basename(f))
                # if os.path.exists(dest_yy) and not overwrite:
                #     pass
                # else:
                #     fix_yy(f, dest)
        else:
            print(f'unknown x, not YY dirs found at: {sys.argv[1]}')


def main():
    if len(sys.argv) < 3:
        print('Usage: python fix_yy.py <src_dir> <dest_dir> [int]')
        sys.exit(1)
    src = sys.argv[1]
    dest = sys.argv[2]
    n_max = int(sys.argv[3]) if len(sys.argv) > 3 else 10 # The latest 10 dirs
    overwrite = False
    fix_yy2(src, dest, n_max, overwrite)


if __name__ == '__main__':
    main()