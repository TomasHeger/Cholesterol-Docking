def classify(pose, group_gap=float('inf')):
    support = float(pose.get('support_fraction', 0))
    inside = float(pose.get('inside_mask_fraction', 0))
    oh = float(pose.get('oh_cosine', 0))
    foreign = float(pose.get('foreign_fraction', 1))
    rescue = bool(pose.get('rescue_pass', False))
    membrane = bool(pose.get('membrane_associated', False))
    ambiguous = group_gap < 0.25

    if (membrane and support >= .93 and inside >= .90 and oh >= .45
            and foreign <= .15 and not ambiguous
            and (not rescue or (support >= .96 and oh >= .70))):
        return 'HIGH', 'strong density, membrane and geometry support'

    if support >= .80 and inside >= .80 and oh >= .15 and foreign <= .30:
        reasons = []
        if rescue: reasons.append('required rescue')
        if ambiguous: reasons.append('neighbor assignment ambiguous')
        if oh < .45: reasons.append('tilted/weak OH direction')
        if foreign > .15: reasons.append('neighboring-density overlap')
        if not membrane: reasons.append('weak membrane association')
        return 'MEDIUM', '; '.join(reasons or ['moderate evidence'])

    reasons = []
    if support < .80: reasons.append('weak density support')
    if inside < .80: reasons.append('partial density-mask occupancy')
    if oh < .15: reasons.append('extreme tilt/weak OH orientation')
    if foreign > .30: reasons.append('substantial neighboring-density overlap')
    if not membrane: reasons.append('not confidently membrane-associated')
    if rescue: reasons.append('required rescue')
    if ambiguous: reasons.append('joint assignment ambiguous')
    return 'LOW', '; '.join(reasons or ['marginal combined evidence'])
