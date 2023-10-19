import base64
import re
import json

class Spreadsheet:
    data: dict[str, str]
    def __init__(self, data):
        self.data = to_dict(data)


    def valid_data(self):
        """
        Get only the non empty rows
        """
        return [row for row in self.data["data"] if len(row) > 0]

    def __str__(self):
        """
        Format data as table in functional style
        """
        return "\n".join(
            [
                "\t".join(self.data["headers"]),
                *["\t".join(row) for row in self.data["data"]]
            ]
        )



    def get_pos(self, key):
        if len(key) != 2:
            raise ValueError("Key should be at most two elements long")
        else:
            x = self.data["headers"].index(key[0])
            y = int(key[1]) - 1
        return y, x

    def __getitem__(self, key):
        x, y = self.get_pos(key)
        return self.data["data"][x][y]

    def __setitem__(self, key, value):
        x, y = self.get_pos(key)
        self.data["data"][x][y] = value


def extract(xml: str) -> str:
    """
    Extract only the "DATA" Portion
    """
    pattern = "<DATA>(.*)<\/DATA>"
    res = re.search(pattern, xml)
    value = res.groups(0)
    return value[0]


def decode_b64(uu_data: bytes) -> str:
    """
    Return the spreadsheet as JSON string
    """
    return base64.b64decode(uu_data)


def to_dict(prop: str) -> dict[str, str]: 
    data = json.loads(decode_b64(extract(prop)))
    return data


def to_spreadsheet(prop: str) -> Spreadsheet: 
    return Spreadsheet(prop)