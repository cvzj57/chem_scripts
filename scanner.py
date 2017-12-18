import subprocess
import sys


class Scanner:
    def __init__(self):
        self.calc_folder_path = '.'
        self.bash_script_name = 'scanup'
        self.distance_range = list(x/10 for x in range(25, 61, 1))
        self.rotation_range = list(x for x in range(0, 360, 10))

    def distance_name_generator(self):
        folder_name_list = []
        for calc in self.distance_range:
            string_parts = str(calc).split('.')
            folder_string = 'd%s_%s' % (string_parts[0], string_parts[1])
            folder_name_list.append(folder_string)
        return folder_name_list

    def rotation_name_generator(self):
        folder_name_list = []
        for calc in self.rotation_range:
            folder_string = 'r%03d' % (calc)
            folder_name_list.append(folder_string)
        return folder_name_list

    def run_scan(self, name_type):
        print('Scan type: ' + str(name_type))
        folder_name_list = None
        if name_type == 'distance':
            print('Scanning range %s to %s...' % (str(self.distance_range[0]), str(self.distance_range[-1])))
            folder_name_list = self.distance_name_generator()
        elif name_type == 'rotation':
            print('Scanning range %s to %s...' % (str(self.rotation_range[0]), str(self.rotation_range[-1])))
            folder_name_list = self.rotation_name_generator()
        for i, folder_name in enumerate(folder_name_list):
            try:
                next_folder_name = folder_name_list[i + 1]
            except IndexError:
                print('End of scan.')
                break
            command = 'bash %s.sh %s %s' % (self.bash_script_name, folder_name, next_folder_name)
            print('running command: %s' % command)
            subprocess.call(command, shell=True, cwd=self.calc_folder_path)


if __name__ == "__main__":
    scanner = Scanner()
    scanner.run_scan(sys.argv[1])
