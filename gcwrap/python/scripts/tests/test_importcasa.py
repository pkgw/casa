import unittest
# Test that casa.py can be imported
class casaimport_test(unittest.TestCase):
    def setUp(self):
        pass
    def test_importcasa(self):
        import casa
    def tearDown(self):
        pass
def suite():
    return [casaimport_test]
