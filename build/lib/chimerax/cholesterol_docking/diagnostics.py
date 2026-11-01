from collections import Counter, defaultdict

class Diagnostics:
    def __init__(self):
        self.sites = defaultdict(lambda: {
            'generated': 0, 'refined': 0, 'rescue_generated': 0,
            'rescue_refined': 0, 'reasons': Counter(), 'best_rejected': None,
            'final_reason': '', 'notes': []
        })

    def generated(self, site, n, rescue=False):
        self.sites[site]['rescue_generated' if rescue else 'generated'] += n

    def refined(self, site, n, rescue=False):
        self.sites[site]['rescue_refined' if rescue else 'refined'] += n

    def reject(self, site, reason, pose=None):
        d = self.sites[site]
        d['reasons'][reason] += 1
        if pose is not None:
            old = d['best_rejected']
            score = float(pose.get('v7_score', pose.get('v6_score', pose.get('score', -1e9))))
            if old is None or score > old[0]:
                d['best_rejected'] = (score, dict(pose))

    def final(self, site, reason):
        self.sites[site]['final_reason'] = reason

    def note(self, site, text):
        self.sites[site]['notes'].append(text)
